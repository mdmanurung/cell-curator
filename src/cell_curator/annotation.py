"""Corpus-driven annotation stages, scheduled by workflow/annotation.smk."""

from __future__ import annotations

import hashlib
import json
import os
import pickle
from html import escape
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .config import CellCuratorConfig, load_config
from .contracts import ContractError, sha256
from .provenance import atomic_publish, publish_json


def annotation_workflow_path() -> Path:
    """Find the same compiled workflow in a checkout or an installed wheel."""
    from importlib.metadata import PackageNotFoundError, distribution

    candidates = [Path(__file__).resolve().parents[2] / "workflow" / "annotation.smk"]
    try:
        dist = distribution("cell-curator")
        candidates.extend(Path(str(dist.locate_file(item))) for item in dist.files or ()
                          if str(item).endswith("share/cell-curator/workflow/annotation.smk"))
    except PackageNotFoundError:
        pass
    snakefile = next((candidate for candidate in candidates if candidate.is_file()), None)
    if snakefile is None:
        raise ContractError("installed corpus annotation Snakefile is missing")
    return snakefile


def execute_annotation_workflow(config_path: str | Path, *, cores: int = 1) -> Path:
    """Execute the compiled DAG through Snakemake, never bypassing its stages."""
    import importlib.util
    import subprocess
    import sys

    from .guided import require_label_gate

    config, root, _ = _setup(config_path)
    require_label_gate(config_path, level=config.guidance.hierarchy_levels[0], parent="ROOT")
    if not isinstance(cores, int) or isinstance(cores, bool) or cores < 1:
        raise ContractError("Snakemake cores must be a positive integer")
    if importlib.util.find_spec("snakemake") is None:
        raise ContractError("the corpus annotation workflow requires the workflow extra (Snakemake)")
    snakefile = annotation_workflow_path()
    command = [sys.executable, "-m", "snakemake", "--snakefile", str(snakefile),
               "--configfile", str(Path(config_path).resolve()), "--cores", str(cores),
               "--directory", str(root)]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode:
        raise ContractError(f"compiled annotation workflow failed:\n{result.stdout}\n{result.stderr}")
    manifest = root / "annotation" / "annotation.manifest.json"
    _check_receipt(manifest)
    return manifest


def annotation_contract(config: CellCuratorConfig) -> dict[str, Any]:
    """Bind rich corpus content and classifier resources before Gate A."""
    from .corpus import load_corpus
    from .knowledge import corpus_context, merge_corpus_context

    development = config.annotation.development_candidates
    if development and not config.annotation.isolation:
        raise ContractError(
            "annotation.development_candidates requires annotation.isolation "
            "(a donor split and a test split), so unreviewed records are checked against the test split"
        )
    corpora = []
    for provider in config.knowledge.providers:
        if provider.kind != "corpus":
            continue
        snapshot = load_corpus(provider.path, expected_sha256=provider.sha256)
        if snapshot.revision != provider.version:
            raise ContractError("configured corpus version differs from its snapshot")
        context = merge_corpus_context(corpus_context(config), provider.context)
        selected = {record.label_id for record in snapshot.select(
            context, development_candidates=development)}
        isolation = None
        if config.annotation.isolation:
            from .corpus import check_test_isolation

            isolation = check_test_isolation(
                snapshot, config.annotation.isolation.donor_split,
                test_split=config.annotation.isolation.test_split,
                dataset_id=config.annotation.isolation.dataset_id)
        corpora.append({"path": provider.path, "sha256": sha256(Path(provider.path)),
                        "isolation": isolation,
                        "revision": snapshot.revision, "context": context,
                        "selected": sorted(selected),
                        "excluded": sorted(record.label_id for record in snapshot.labels
                                           if record.label_id not in selected)})
    resources = {}
    for classifier in config.annotation.classifiers:
        files: dict[str, str] = {}
        for kind in ("model", "reference", "python_executable"):
            value = getattr(classifier, kind)
            if not value:
                continue
            path = Path(value)
            if path.is_file():
                files[str(path)] = sha256(path)
            elif path.is_dir():
                files.update({str(item): sha256(item) for item in sorted(path.rglob("*"))
                              if item.is_file()})
            elif classifier.required:
                raise ContractError(f"required classifier {kind} is absent: {path}")
        resources[classifier.name] = files
    return {"corpora": corpora, "classifier_resources": resources,
            "development_candidates": development,
            "rules": config.annotation.model_dump(mode="json")}


def _setup(config_path: str | Path) -> tuple[CellCuratorConfig, Path, dict[str, Any]]:
    from .guided import load_state
    from .pipeline import _programs, initialize, run_root

    config = load_config(config_path)
    if not config.annotation.enabled:
        raise ContractError("the corpus annotation DAG requires annotation.enabled: true")
    root = run_root(config)
    state = load_state(root)
    if state.config_sha256 != sha256(Path(config_path)):
        raise ContractError("annotation configuration changed; use a new run_id")
    initialize(Path(config_path))
    plan = json.loads((root / "plan.json").read_text())
    if plan.get("annotation_contract") != annotation_contract(config):
        raise ContractError("corpus, classifier resources, or rules changed after planning")
    _programs(config)
    return config, root, plan


def _input(config: CellCuratorConfig, root: Path) -> Any:
    from .data import load_input

    manifest = json.loads((root / "input_contract.manifest.json").read_text())
    frozen = config.model_copy(update={"input": config.input.model_copy(
        update={"path": manifest["canonical_path"]})})
    return load_input(frozen)



def _ranking_positions(groups: pd.Series, cap: int | None, seed: int) -> np.ndarray:
    """Positions of the cells the de novo ranking sees: at most `cap` per cluster.

    Cells are sorted by name before each cluster's draw, and each cluster has its own
    generator, so a cluster's subset does not depend on the other clusters in the call.
    """

    if cap is None:
        return np.arange(len(groups))
    names, labels = groups.index.astype(str).to_numpy(), groups.astype(str).to_numpy()
    keep = []
    for cluster in sorted(set(labels)):
        members = np.flatnonzero(labels == cluster)
        members = members[np.argsort(names[members], kind="stable")]
        if len(members) > cap:
            seed_words = [seed, int(hashlib.sha256(cluster.encode()).hexdigest()[:8], 16)]
            members = members[np.sort(np.random.default_rng(seed_words).choice(len(members), cap, replace=False))]
        keep.append(members)
    return np.sort(np.concatenate(keep))


def _ranked_markers(matrix: Any, groups: pd.Series, features: pd.Index, detection: Any, config: CellCuratorConfig,
                    context: dict[str, Any], role: str) -> pd.DataFrame:
    """compute_bottom_up_markers on the capped cells, read from or written to the ranking cache."""

    from . import __version__
    from . import evidence as evidence_module
    from .evidence import compute_bottom_up_markers

    settings = config.markers
    marker_config = settings.model_dump(mode="json")
    positions = _ranking_positions(groups, settings.ranking_max_cells_per_cluster, settings.ranking_seed)
    ranked = groups.iloc[positions]
    counts = {str(k): int(v) for k, v in ranked.astype(str).value_counts().sort_index().items()}
    key = hashlib.sha256(json.dumps({
        "input_sha256": context["input_sha256"], "role": role, "version": __version__,
        "ranking_code_sha256": sha256(Path(evidence_module.__file__)), "marker_config": marker_config,
        "detection": context["detection"], "features": [str(item) for item in features],
        "cells": [[str(name), str(cluster)] for name, cluster in zip(ranked.index, ranked, strict=True)],
    }, sort_keys=True).encode()).hexdigest()
    cache = Path(settings.ranking_cache_dir) / f"{key}.pkl" if settings.ranking_cache_dir else None
    frame, status = None, "off" if cache is None else "miss"
    if cache is not None and cache.is_file():
        with cache.open("rb") as handle:
            stored = pickle.load(handle)
        if stored["key"] == key:
            frame, status = stored["frame"], "hit"
    if frame is None:
        frame = compute_bottom_up_markers(matrix[positions], ranked, features, top_n=None, methods=("wilcoxon", "logreg"),
            detection=detection[positions], marker_config=marker_config)
        if cache is not None:
            cache.parent.mkdir(parents=True, exist_ok=True)
            partial = cache.with_name(f".{cache.name}.partial.{os.getpid()}")
            partial.write_bytes(pickle.dumps({"key": key, "frame": frame}, protocol=4))
            os.replace(partial, cache)
    frame.attrs["runtime"] = {**frame.attrs.get("runtime", {}), "ranking_cache": status, "ranking_cache_key": key,
        "ranking_max_cells_per_cluster": settings.ranking_max_cells_per_cluster,
        "ranking_seed": settings.ranking_seed, "ranked_cells_per_cluster": counts}
    return frame


def _frame(path: Path, frame: pd.DataFrame) -> None:
    atomic_publish(path, frame.to_csv(sep="\t", index=False, lineterminator="\n").encode())


def _receipt(path: Path, *, artifacts: list[Path], inputs: list[Path], **values: Any) -> dict[str, Any]:
    value = {"schema_version": 3, **values,
             "inputs": {str(item): sha256(item) for item in inputs},
             "artifacts": {str(item): sha256(item) for item in artifacts}}
    publish_json(path, value)
    return value


def _check_receipt(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ContractError(f"annotation receipt is missing: {path}")
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ContractError(f"invalid annotation receipt: {path}") from exc
    if not isinstance(value, dict) or any(not isinstance(value.get(key), dict) for key in ("inputs", "artifacts")):
        raise ContractError(f"annotation receipt lacks input/artifact hashes: {path}")
    for key in ("inputs", "artifacts"):
        for name, digest in value[key].items():
            if not Path(name).is_file() or sha256(Path(name)) != digest:
                raise ContractError(f"stale annotation artifact: {name}")
    return value


def run_classifier(config_path: str | Path, name: str) -> dict[str, Any]:
    from .crosschecks import celltypist_hypothesis, popv_hypothesis
    from .data import close_input
    from .guided import require_label_gate

    config, root, _ = _setup(config_path)
    require_label_gate(config_path, level=config.guidance.hierarchy_levels[0], parent="ROOT")
    classifier = next((item for item in config.annotation.classifiers if item.name == name), None)
    if classifier is None:
        raise ContractError(f"classifier is not configured: {name}")
    target = root / "annotation" / "classifiers" / name
    if classifier.python_executable:
        import subprocess
        import sys

        executable = Path(classifier.python_executable)
        # Distinct virtualenv executables often symlink to the same Python binary.
        if os.path.abspath(executable) != os.path.abspath(sys.executable):
            if not executable.is_file():
                raise ContractError(f"classifier Python interpreter is absent: {executable}")
            result = subprocess.run([str(executable), "-m", "cell_curator", "annotation", "classifier",
                "--config", str(Path(config_path).resolve()), "--name", name],
                text=True, capture_output=True, check=False)
            if result.returncode:
                raise ContractError(f"isolated classifier {name} failed:\n{result.stdout}\n{result.stderr}")
            status_path = target / "status.json"
            if not status_path.is_file():
                raise ContractError("isolated classifier did not publish its status receipt")
            try:
                status = json.loads(status_path.read_text())
            except (OSError, ValueError) as exc:
                raise ContractError("isolated classifier published an invalid status receipt") from exc
            if not isinstance(status, dict) or status.get("classifier") != name:
                raise ContractError("isolated classifier status names a different classifier")
            cells_path = target / "cell_predictions.tsv"
            if not cells_path.is_file() or sha256(cells_path) != status.get("prediction_sha256"):
                raise ContractError("isolated classifier evidence is absent or changed")
            if classifier.required and status.get("status") != "complete":
                raise ContractError(f"required isolated classifier {name} did not complete")
            return status
    loaded = _input(config, root)
    try:
        adata = loaded.assays[config.annotation.assay].copy()
        cluster_key = config.input.keys.canonical_parent
        adata.obs[cluster_key] = loaded.root.obs.loc[adata.obs_names, cluster_key].astype(str)
        if classifier.method == "celltypist":
            if classifier.counts_layer not in adata.layers:
                raise ContractError("CellTypist requires the declared counts layer")
            counts = adata.layers[classifier.counts_layer].copy()
            values = counts.data if sparse.issparse(counts) else np.asarray(counts)
            if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.round(values)):
                raise ContractError("CellTypist counts must be finite, nonnegative integers")
            totals = np.asarray(counts.sum(axis=1)).ravel()
            scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
            normalized = sparse.diags(scale) @ counts if sparse.issparse(counts) else counts * scale[:, None]
            if sparse.issparse(normalized):
                normalized = normalized.tocsr()
                normalized.data = np.log1p(normalized.data)
            else:
                normalized = np.log1p(normalized)
            adata.X = normalized
            cells, status = celltypist_hypothesis(adata, model=classifier.model,
                cluster_column=cluster_key, output_dir=target)
        else:
            import anndata as ad

            reference = ad.read_h5ad(classifier.reference)
            cells, status = popv_hypothesis(adata, reference, model_dir=classifier.model,
                reference_label=classifier.reference_label, methods=tuple(classifier.methods),
                query_counts_layer=classifier.counts_layer,
                reference_counts_layer=classifier.reference_counts_layer,
                cluster_column=cluster_key, output_dir=target, seed=config.annotation.seed)
        if cells.empty:
            cells = pd.DataFrame(columns=["barcode", "predicted_label", "cluster_id"])
        _frame(target / "cell_predictions.tsv", cells)
        status = {**status, "required": classifier.required, "classifier": name,
                  "prediction_sha256": sha256(target / "cell_predictions.tsv")}
        publish_json(target / "status.json", status)
    finally:
        close_input(loaded)
    if classifier.required and status["status"] != "complete":
        raise ContractError(f"required classifier {name} did not complete: {status.get('reason', '')}")
    return status


def _predictions(
    config: CellCuratorConfig, root: Path
) -> tuple[dict[str, pd.DataFrame], list[Path], dict[str, set[str]]]:
    predictions, paths = {}, []
    dependent: dict[str, set[str]] = {}
    for classifier in config.annotation.classifiers:
        directory = root / "annotation" / "classifiers" / classifier.name
        status_path, cells_path = directory / "status.json", directory / "cell_predictions.tsv"
        if not status_path.is_file() or not cells_path.is_file():
            raise ContractError(f"classifier evidence is missing: {classifier.name}; run its Snakemake stage")
        status = json.loads(status_path.read_text())
        if sha256(cells_path) != status["prediction_sha256"]:
            raise ContractError("classifier predictions changed after publication")
        paths.extend([status_path, cells_path])
        if status["status"] != "complete":
            if classifier.required:
                raise ContractError(f"required classifier is unavailable: {classifier.name}")
            continue
        predictions[classifier.name] = pd.read_csv(cells_path, sep="\t", keep_default_na=False, dtype={"cluster_id": str, "barcode": str})
        dependent[classifier.name] = {str(item).casefold() for item in status.get("dependent_methods", [])}
    return predictions, paths, dependent


def _signature_support(
    detection: Any, groups: pd.Series, feature_lookup: dict[str, int],
    positive_genes: tuple[str, ...], negative_genes: tuple[str, ...], config: CellCuratorConfig,
) -> tuple[dict[str, float], bool, float]:
    """Decision row 17 for one gene signature: per cluster, the share of cells that detect
    enough positive genes and no negative gene. Returns (support by cluster, measured, coverage)."""
    pos = [feature_lookup[gene.casefold()] for gene in positive_genes if gene.casefold() in feature_lookup]
    neg = [feature_lookup[gene.casefold()] for gene in negative_genes if gene.casefold() in feature_lookup]
    coverage = len(pos) / max(len(positive_genes), 1)
    negative_coverage = len(neg) / len(negative_genes) if negative_genes else 1.
    minimum = config.evidence.min_positive_coverage
    measured = bool(pos) and coverage >= minimum and negative_coverage >= minimum
    support: dict[str, float] = {}
    for cluster in sorted(set(groups)):
        if not measured:
            support[cluster] = float("nan")
            continue
        values = detection[groups.eq(cluster).to_numpy()]
        positive = values[:, pos].toarray() if sparse.issparse(values) else np.asarray(values)[:, pos]
        negative = values[:, neg].toarray() if sparse.issparse(values) else np.asarray(values)[:, neg]
        cell_positive = (positive > 0).mean(axis=1) >= config.annotation.min_state_positive_fraction
        support[cluster] = float((cell_positive & ~(negative > 0).any(axis=1)).mean())
    return support, measured, coverage


def _program_exclusions(config: CellCuratorConfig, programs: list[Any]) -> dict[str, tuple[Any, ...]]:
    """Decision row 30: the exclusions each corpus-backed program carries, by program source."""
    from .corpus import load_corpus

    found: dict[str, tuple[Any, ...]] = {}
    for provider in config.knowledge.providers:
        if provider.kind != "corpus":
            continue
        snapshot = load_corpus(provider.path, expected_sha256=provider.sha256)
        by_id = {label.label_id: label for label in snapshot.labels}
        prefix = f"corpus:{snapshot.corpus_id}:{snapshot.digest}:"
        for program in programs:
            if program.source.startswith(prefix) and program.source[len(prefix):] in by_id:
                found[program.source] = by_id[program.source[len(prefix):]].exclusions
    return found


def _marker_family(frames: list[pd.DataFrame], adjusted_pvalue_max: float, top_n: int = 100) -> pd.DataFrame:
    """Decision rows 7 and 8: one BH family per parent scope and method, then the top
    genes by effect among the tests that pass it."""
    from .evidence import benjamini_hochberg

    frame = pd.concat(frames, ignore_index=True)
    wilcoxon = frame["method"].eq("wilcoxon")
    frame.loc[wilcoxon, "adjusted_pvalue"] = benjamini_hochberg(frame.loc[wilcoxon, "pvalue"].to_numpy())
    frame = frame[~wilcoxon | frame["adjusted_pvalue"].le(adjusted_pvalue_max)]
    keys = ["cluster_id", "comparator_role", "method"]
    frame = frame.sort_values(keys + ["effect", "feature"], ascending=[True, True, True, False, True])
    frame = frame.groupby(keys, sort=False).head(top_n).copy()
    frame["rank"] = frame.groupby(keys, sort=False).cumcount() + 1
    return frame.reset_index(drop=True)


def _is_coarser(label: str, candidate: str, parents: dict[str, str]) -> bool:
    """True when ``label`` is a strict ancestor of ``candidate``."""
    while candidate in parents:
        candidate = parents[candidate]
        if candidate == label:
            return True
    return False


def _classifier_counts(raw: pd.Series, label_map: dict[str, str], canonical: set[str],
                       parents: dict[str, str], candidate: str) -> dict[str, int]:
    """Decision rows 12 and 13: a coarser label abstains, an unmapped label is dropped."""
    counts = {"agree": 0, "disagree": 0, "coarser": 0, "unmapped": 0}
    for value in raw:
        label = value if value in canonical else label_map.get(value)
        if label is None:
            counts["unmapped"] += 1
        elif _compatible(label, candidate, parents):
            counts["agree"] += 1
        elif _is_coarser(label, candidate, parents):
            counts["coarser"] += 1
        else:
            counts["disagree"] += 1
    return counts


def _compatible(label: str, candidate: str, parents: dict[str, str]) -> bool:
    while label in parents:
        if label == candidate:
            return True
        label = parents[label]
    return label == candidate


def collect_level(config_path: str | Path, level: str) -> dict[str, Any]:
    from .capabilities import build_capability_registry
    from .data import close_input, feature_ids_for, matrix_for
    from .evidence import score_marker_programs
    from .guided import require_label_gate
    from .markers import marker_dotplot_svg, profile_marker_panel
    from .pipeline import _programs

    config, root, plan = _setup(config_path)
    levels = config.guidance.hierarchy_levels
    if level not in levels:
        raise ContractError(f"unknown annotation level: {level}")
    position = levels.index(level)
    target = root / "annotation" / "levels" / level
    vocabulary = pd.read_csv(plan["vocabulary_path"], sep="\t", keep_default_na=False, dtype={"cluster_id": str, "barcode": str})
    identities = vocabulary[vocabulary["axis"].eq("cell_type")]
    scoped = identities[identities["level"].eq(level)]
    previous: dict[str, str] = {}
    inputs = [Path(config_path), root / "plan.json", Path(plan["vocabulary_path"])]
    if position:
        previous_dir = root / "annotation" / "levels" / levels[position - 1]
        _check_receipt(previous_dir / "decision.manifest.json")
        prior = pd.read_csv(previous_dir / "labels.tsv", sep="\t", keep_default_na=False, dtype={"cluster_id": str, "barcode": str})
        previous = dict(zip(prior["cluster_id"].astype(str), prior["label"].astype(str), strict=True))
        inputs.append(previous_dir / "decision.manifest.json")
    predictions, classifier_paths, dependent = _predictions(config, root)
    inputs.extend(classifier_paths)
    loaded = _input(config, root)
    try:
        adata = loaded.assays[config.annotation.assay]
        assay = config.input.assays[config.annotation.assay]
        declaration = assay.representations[config.annotation.representation]
        if declaration.ranking_backend.lower() not in {"scipy_mannwhitney", "scipy", "cpu"}:
            raise ContractError("the corpus DAG currently requires an explicit CPU marker backend; native GPU lanes remain in workflow/Snakefile")
        matrix = matrix_for(adata, declaration)
        detection = matrix_for(adata, declaration.detection.model_dump()) if declaration.detection else None
        if detection is None or detection.shape != matrix.shape:
            raise ContractError("corpus annotation requires aligned detection and ranking matrices")
        features = feature_ids_for(adata, assay.feature_id_column, declaration=declaration,
                                   gene_mapping_path=assay.gene_mapping_path)
        groups = loaded.root.obs.loc[adata.obs_names, config.input.keys.canonical_parent].astype(str)
        groups.index = adata.obs_names.astype(str)
        capability = next(item for item in build_capability_registry(config)
                          if item.assay == config.annotation.assay and item.representation == config.annotation.representation)
        detection_features = feature_ids_for(adata, assay.feature_id_column, declaration=declaration.detection,
                                             gene_mapping_path=assay.gene_mapping_path)
        if not features.equals(detection_features):
            raise ContractError("detection and ranking feature orders differ")
        programs = [item for item in _programs(config) if item.feature_space == capability.feature_space]
        feature_lookup = {str(item).casefold(): index for index, item in enumerate(features)}
        exclusions_by_source = _program_exclusions(config, programs)
        exclusion_cache: dict[tuple[Any, ...], tuple[dict[str, float], bool, float]] = {}
        score_frames: list[pd.DataFrame] = []
        marker_frames: list[pd.DataFrame] = []
        statuses: list[dict[str, Any]] = []
        marker_runtime = []
        ranking_context = {"input_sha256": sha256(Path(json.loads((root / "input_contract.manifest.json").read_text())["canonical_path"])),
                           "detection": declaration.detection.model_dump(mode="json") if declaration.detection else None}
        for parent, block in scoped.groupby("parent_scope", sort=True):
            parent = str(parent)
            allowed = set(groups) if position == 0 else {cluster for cluster, label in previous.items() if label == parent}
            mask = groups.isin(allowed).to_numpy()
            if not mask.any():
                statuses.append({"parent": parent, "status": "unreachable"})
                continue
            require_label_gate(config_path, level=level, parent=parent)
            scoped_groups = groups.iloc[np.flatnonzero(mask)]
            candidate_programs = [item for item in programs if item.axis == "cell_type"
                                  and item.label in set(block["label"]) and (item.parent or "ROOT") == parent]
            scores = score_marker_programs(matrix[mask], scoped_groups, features, candidate_programs,
                detection=detection[mask], run_id=config.run.run_id, assay=config.annotation.assay,
                representation=config.annotation.representation, capability=capability,
                min_negative_coverage=config.evidence.min_positive_coverage)
            scored_names = set(zip(scores["cluster_id"], scores["label"], strict=True))
            absent = []
            for cluster in sorted(set(scoped_groups)):
                for program in candidate_programs:
                    if (cluster, program.label) not in scored_names:
                        absent.append({"cluster_id": cluster, "label": program.label, "score": 0.,
                            "vote_fraction": 0., "program_coverage": 0., "positive_coverage": 0.,
                            "negative_violation": float("nan"), "detection_valid": True})
            if absent:
                scores = pd.concat([scores, pd.DataFrame(absent)], ignore_index=True)
            clusters = sorted(set(scoped_groups))
            bottom = pd.DataFrame()
            first_frame = len(marker_frames)
            comparators_ran: set[tuple[str, str]] = set()
            if len(clusters) > 1:
                bottom = _ranked_markers(matrix[mask], scoped_groups, features, detection[mask], config,
                                         ranking_context, "parent_remainder")
                comparators_ran.update((cluster, "parent_remainder") for cluster in clusters)
                marker_runtime.append({"parent": parent, "role": "parent_remainder", **bottom.attrs.get("runtime", {})})
                bottom["comparator_role"] = "parent_remainder"
                marker_frames.append(bottom.assign(level=level, parent_scope=parent))
                embedding = np.asarray(loaded.root.obsm[config.input.keys.embedding])[mask]
                centroids = {cluster: embedding[scoped_groups.to_numpy() == cluster].mean(axis=0) for cluster in clusters}
                provisional = scores.sort_values(["vote_fraction", "score", "label"], ascending=[False, False, True]).drop_duplicates("cluster_id").set_index("cluster_id")["label"].to_dict()
                for cluster in clusters:
                    others = [item for item in clusters if item != cluster]
                    # Decision row 5: prefer a cluster with a different provisional label;
                    # when every cluster shares one, fall back to the nearest cluster.
                    competitors = [item for item in others if provisional[item] != provisional[cluster]]
                    selection = "nearest_cluster_with_different_provisional_marker_identity"
                    if not competitors:
                        competitors, selection = others, "nearest_cluster_all_share_provisional_identity"
                    sibling = min(competitors,
                                  key=lambda item: (float(np.linalg.norm(centroids[cluster] - centroids[item])), item))
                    pair = scoped_groups.isin([cluster, sibling]).to_numpy()
                    pair_markers = _ranked_markers(matrix[mask][pair], scoped_groups.iloc[np.flatnonzero(pair)],
                        features, detection[mask][pair], config, ranking_context, "closest_sibling")
                    comparators_ran.add((cluster, "closest_sibling"))
                    marker_runtime.append({"parent": parent, "cluster": cluster, "role": "closest_sibling",
                                           **pair_markers.attrs.get("runtime", {})})
                    pair_markers = pair_markers[pair_markers["cluster_id"].eq(cluster)].copy()
                    pair_markers["comparator_role"] = "closest_sibling"
                    pair_markers["reference_cluster"] = sibling
                    pair_markers["comparator_selection"] = selection
                    marker_frames.append(pair_markers.assign(level=level, parent_scope=parent))
                statuses.append({"parent": parent, "status": "complete", "marker_methods": ["wilcoxon", "logreg"]})
            else:
                statuses.append({"parent": parent, "status": "unavailable_no_comparator", "marker_methods": []})
            if len(marker_frames) > first_frame:
                marker_frames[first_frame:] = [_marker_family(
                    marker_frames[first_frame:], config.evidence.adjusted_pvalue_max)]
            scope_markers = marker_frames[first_frame] if len(marker_frames) > first_frame else pd.DataFrame()
            scores["raw_vote_fraction"] = scores["vote_fraction"]
            for index, row in scores.iterrows():
                program = next(item for item in candidate_programs if item.label == row["label"])
                failures = []
                for field in ("program_coverage", "positive_coverage"):
                    if not np.isfinite(row[field]) or row[field] < config.evidence.min_positive_coverage:
                        failures.append(field)
                # Decision row 2: unmeasured negatives are not contradicted; flag, do not fail.
                flags = []
                if not np.isfinite(row["negative_violation"]):
                    flags.append("unmeasured_negatives")
                elif row["negative_violation"] > config.evidence.max_negative_violation:
                    failures.append("negative_markers")
                if not program.negative:
                    flags.append("no_declared_negatives")  # decision row 4
                for excluded in exclusions_by_source.get(program.source, ()):
                    # Decision row 30: a population the label must not match fails it when the
                    # cluster shows that population's own signature (the row 17 state rule).
                    signature = (excluded.name, excluded.positive, excluded.negative)
                    if signature not in exclusion_cache:
                        exclusion_cache[signature] = _signature_support(
                            detection, groups, feature_lookup, excluded.positive, excluded.negative, config)
                    by_cluster, excluded_measured, _ = exclusion_cache[signature]
                    if not excluded_measured:
                        flags.append(f"exclusion_unmeasured:{excluded.name}")
                    elif by_cluster[str(row["cluster_id"])] >= config.annotation.min_state_fraction:
                        failures.append(f"excluded:{excluded.name}")
                # Decision row 6: at least 2 positive markers, or all if the label has 1.
                needed = min(2, len(program.positive))
                for role in ("parent_remainder", "closest_sibling"):
                    if (str(row["cluster_id"]), role) not in comparators_ran:
                        flags.append(f"comparator_unavailable:{role}")  # decision row 3
                        continue
                    for method in ("wilcoxon", "logreg"):
                        evidence = scope_markers[
                            scope_markers["cluster_id"].eq(row["cluster_id"])
                            & scope_markers["method"].eq(method) & scope_markers["comparator_role"].eq(role)
                            & scope_markers["effect"].gt(0)]
                        hits = set(evidence["feature"].astype(str).str.casefold()).intersection(
                            gene.casefold() for gene in program.positive)
                        if len(hits) < needed:
                            failures.append(f"{role}:{method}")
                classifier_counts: dict[str, dict[str, int]] = {}
                agreement_by_name: dict[str, float] = {}
                parents = dict(zip(identities["label"], identities["parent_scope"], strict=True))
                for classifier in config.annotation.classifiers:
                    if classifier.name not in predictions:
                        continue
                    cells = predictions[classifier.name]
                    expected = set(groups.index[groups.eq(str(row["cluster_id"]))])
                    if cells["barcode"].duplicated().any() or set(cells["barcode"].astype(str)) != set(groups.index):
                        raise ContractError("classifier barcodes differ from the annotation input")
                    raw = cells[cells["barcode"].astype(str).isin(expected)]["predicted_label"].astype(str)
                    counts = _classifier_counts(raw, classifier.label_map, set(identities["label"]), parents, program.label)
                    classifier_counts[classifier.name] = counts
                    judged = counts["agree"] + counts["disagree"]
                    agreement_by_name[classifier.name] = counts["agree"] / judged if judged else float("nan")
                # Decision row 11: a classifier that runs another one inside it counts as one vote with it.
                classifier_support: dict[str, float] = {}
                for classifier in config.annotation.classifiers:
                    if classifier.name not in agreement_by_name:
                        continue
                    inside = dependent.get(classifier.name, set())
                    members = sorted({classifier.name} | {other.name for other in config.annotation.classifiers
                                                          if other.name in agreement_by_name and other.method in inside}
                                     | {other.name for other in config.annotation.classifiers
                                        if other.name in agreement_by_name and classifier.method in dependent.get(other.name, set())})
                    key = "+".join(members)
                    if key in classifier_support:
                        continue
                    values = [agreement_by_name[name] for name in members if np.isfinite(agreement_by_name[name])]
                    if not values:
                        continue
                    classifier_support[key] = float(np.mean(values))
                    if classifier_support[key] < config.annotation.min_classifier_agreement:
                        # Decision row 10: advisory by default; `required` makes it a gate again.
                        if any(item.required for item in config.annotation.classifiers if item.name in members):
                            failures.append(f"classifier:{key}")
                        else:
                            flags.append(f"classifier_disagreement:{key}")
                scores.at[index, "gate_failures"] = "|".join(failures)
                scores.at[index, "gate_flags"] = "|".join(flags)
                scores.at[index, "classifier_support"] = json.dumps(classifier_support, sort_keys=True)
                scores.at[index, "classifier_counts"] = json.dumps(classifier_counts, sort_keys=True)
                if failures:
                    scores.at[index, "vote_fraction"] = 0.
            scores["level"], scores["parent_scope"] = level, parent
            score_frames.append(scores)
        scores = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame(columns=["cluster_id", "label", "score", "vote_fraction", "level", "parent_scope"])
        markers = pd.concat(marker_frames, ignore_index=True) if marker_frames else pd.DataFrame(columns=["cluster_id", "feature", "method", "comparator_role", "level", "parent_scope"])
        # States are tested independently, so co-occurring programs can both pass.
        state_rows = []
        for program in (item for item in programs if item.axis == "cell_state"):
            support_by_cluster, measured, coverage = _signature_support(
                detection, groups, feature_lookup, program.positive, program.negative, config)
            for cluster, support in support_by_cluster.items():
                state_rows.append({"cluster_id": cluster, "state": program.label, "parent_scope": program.parent or "ROOT",
                    "fraction": support, "measured": measured, "coverage": coverage,
                    "passes": measured and support >= config.annotation.min_state_fraction})
        states = pd.DataFrame(state_rows, columns=["cluster_id", "state", "parent_scope", "fraction", "measured", "coverage", "passes"])
        profile = profile_marker_panel(matrix, groups, features, programs, detection=detection)
        plot = marker_dotplot_svg(profile)
        for name, frame in (("scores.tsv", scores), ("markers.tsv", markers), ("states.tsv", states)):
            _frame(target / name, frame)
        atomic_publish(target / "marker_evidence.svg", plot.encode())
    finally:
        close_input(loaded)
    return _receipt(target / "evidence.manifest.json", artifacts=[target / name for name in ("scores.tsv", "markers.tsv", "states.tsv", "marker_evidence.svg")],
                    inputs=inputs, artifact_kind="cell_curator_corpus_level_evidence", level=level,
                    scopes=statuses, marker_runtime=marker_runtime,
                    ranking={"max_cells_per_cluster": config.markers.ranking_max_cells_per_cluster,
                             "seed": config.markers.ranking_seed, "cache_dir": config.markers.ranking_cache_dir},
                    unit_of_inference="descriptive_cluster_markers")


def finalize_annotation(config_path: str | Path) -> dict[str, Any]:
    from .data import close_input
    from .guided import require_label_gate
    from .pipeline import execute_recursive_hierarchy

    config, root, plan = _setup(config_path)
    score_tables, decision_tables, inputs = {}, {}, []
    decisions = []
    for level in config.guidance.hierarchy_levels:
        directory = root / "annotation" / "levels" / level
        _check_receipt(directory / "evidence.manifest.json")
        _check_receipt(directory / "decision.manifest.json")
        table = pd.read_csv(directory / "scores.tsv", sep="\t", keep_default_na=False, dtype={"cluster_id": str, "barcode": str})
        for parent, block in table.groupby("parent_scope", sort=True):
            score_tables[level, str(parent)] = block
        labels = pd.read_csv(directory / "labels.tsv", sep="\t", keep_default_na=False, dtype={"cluster_id": str, "barcode": str})
        for parent, block in labels.groupby("parent_scope", sort=True):
            decision_tables[level, str(parent)] = block
        decisions.append(labels)
        inputs.extend([directory / "decision.manifest.json", directory / "evidence.manifest.json"])
    completed = root / "annotation" / "annotation.manifest.json"
    if completed.is_file():
        return _check_receipt(completed)
    execution = execute_recursive_hierarchy(config_path, score_tables=score_tables,
                                            decisions=decision_tables)
    loaded = _input(config, root)
    try:
        cells = pd.DataFrame({"barcode": loaded.root.obs_names.astype(str),
            "cluster_id": loaded.root.obs[config.input.keys.canonical_parent].astype(str).to_numpy()})
    finally:
        close_input(loaded)
    leaf = pd.Series("Unknown", index=cells.index)
    ancestry: dict[str, set[str]] = {cluster: set() for cluster in cells["cluster_id"].unique()}
    for level, table in zip(config.guidance.hierarchy_levels, decisions, strict=True):
        mapping = dict(zip(table["cluster_id"].astype(str), table["label"].astype(str), strict=True))
        current = cells["cluster_id"].map(mapping).fillna("Unknown")
        leaf = current.where(current.ne("Unknown"), leaf)
        cells[level] = leaf
        for cluster, label in mapping.items():
            if label != "Unknown":
                ancestry[cluster].add(label)
    states = pd.read_csv(root / "annotation" / "levels" / config.guidance.hierarchy_levels[-1] / "states.tsv", sep="\t", dtype={"cluster_id": str})
    vocabulary = pd.read_csv(plan["vocabulary_path"], sep="\t", keep_default_na=False)
    state_vocabulary = vocabulary[vocabulary["axis"].eq("cell_state")]
    state_map = {}
    for cluster, block in states.groupby("cluster_id", sort=True):
        accepted = block[block["passes"].astype(str).str.lower().eq("true")]
        selected = set()
        for row in accepted.itertuples():
            parent = str(row.parent_scope)
            if parent != "ROOT" and parent not in ancestry.get(str(cluster), set()):
                continue
            permitted = state_vocabulary[state_vocabulary["label"].eq(str(row.state))
                                         & state_vocabulary["parent_scope"].eq(parent)]
            if len(permitted) != 1:
                raise ContractError("state evidence is absent or ambiguous in the approved vocabulary")
            require_label_gate(config_path, level=str(permitted.iloc[0]["level"]), parent=parent)
            selected.add(str(row.state))
        state_map[str(cluster)] = sorted(selected)
    cells["cell_state"] = cells["cluster_id"].map(lambda cluster: "|".join(state_map.get(cluster, [])))
    cells["L4"] = leaf + cells["cell_state"].map(lambda value: f" {value.replace('|', ' ')}" if value else "")
    cells["annotation_unit"] = "cluster_propagated_to_cells"
    development = bool(plan["annotation_contract"].get("development_candidates"))
    if development:
        cells["development_candidates"] = True
    target = root / "annotation"
    _frame(target / "cell_labels.tsv", cells)
    mapping = pd.DataFrame({"barcode": cells["barcode"], "canonical_parent": cells["cluster_id"],
                            "leaf_id": cells["cluster_id"]})
    if mapping["barcode"].duplicated().any() or len(mapping) != len(cells):
        raise ContractError("corpus mapping must contain exactly one row per source cell")
    mapping_path = root / "mapping" / "parent_child_mapping.tsv"
    _frame(mapping_path, mapping)
    publish_json(root / "mapping" / "mapping.validation.json", {
        "schema_version": 3, "artifact_kind": "cell_curator_final_mapping_validation",
        "status": "complete", "mapping_sha256": sha256(mapping_path),
        "n_cells": len(mapping), "n_leaves": int(mapping["leaf_id"].nunique()),
        "invariants": {"one_row_per_input_cell": not mapping["barcode"].duplicated().any(),
                       "canonical_parent_membership_unchanged": mapping["canonical_parent"].equals(cells["cluster_id"]),
                       "parent_counts_conserved": mapping["canonical_parent"].value_counts().equals(cells["cluster_id"].value_counts())}})
    report = ["<!doctype html><html lang='en'><meta charset='utf-8'><title>Annotation evidence</title>",
        *(["<p><strong>DEVELOPMENT CANDIDATES: unreviewed corpus records were called. "
           "This run can never be written back.</strong></p>"] if development else []),
        "<h1>Corpus annotation evidence</h1><p>Cluster labels propagated to cells. Classifiers provide corroboration; confidence is not a calibrated probability.</p>"]
    for level, table in zip(config.guidance.hierarchy_levels, decisions, strict=True):
        report.extend([f"<h2>{escape(level)}</h2>", table.to_html(index=False, escape=True),
                       f"<img alt='Measured marker expression for {escape(level)}' src='levels/{escape(level)}/marker_evidence.svg'>",
                       f"<p><a href='levels/{escape(level)}/scores.tsv'>Candidate evidence and failed gates</a></p>"])
    report.append("<p>Human assumptions, configured critic passes, and write-back approval remain required.</p></html>")
    atomic_publish(target / "report.html", "\n".join(report).encode())
    return _receipt(target / "annotation.manifest.json", artifacts=[target / "cell_labels.tsv", target / "report.html",
        mapping_path, root / "mapping" / "mapping.validation.json",
        *sorted((root / "labels").glob("L*.tsv"))],
        inputs=[*inputs, root / "hierarchy" / "execution.manifest.json"],
        artifact_kind="cell_curator_corpus_annotation", status="complete", corpus=plan["annotation_contract"]["corpora"],
        hierarchy_validation=execution["hierarchy_validation"],
        annotation_unit="cluster_propagated_to_cells", writeback_applied=False,
        development_candidates=development)
