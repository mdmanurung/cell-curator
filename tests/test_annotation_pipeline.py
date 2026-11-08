from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import REPO_ROOT
from test_end_to_end import _config

from cell_curator.annotation import (
    annotation_contract,
    collect_level,
    finalize_annotation,
    resolve_level,
    run_classifier,
)
from cell_curator.config import load_config
from cell_curator.contracts import ContractError, sha256
from cell_curator.corpus import CorpusLabel, CorpusSnapshot, CorpusSource
from cell_curator.guided import confirm_context, construct_plan, prepare_scene
from cell_curator.provenance import publish_json


def annotation_fixture(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    classifier: bool = False,
    conflict: bool = False,
    unavailable: bool = False,
    approved: bool = True,
    collapse_lineage: bool = False,
    required: bool = True,
    exclusions: dict[str, list[dict]] | None = None,
) -> tuple[Path, Path, dict[str, list[str]]]:
    """Eight terminal clusters share their broad markers with sibling clusters."""
    features = [
        "lineage-A",
        "lineage-B",
        *[f"branch-{i}" for i in range(4)],
        *[f"leaf-{i}" for i in range(8)],
        "state-S",
    ]
    counts = np.zeros((160, len(features)), dtype=np.float32)
    expected: dict[str, list[str]] = {}
    labels = []

    def label(
        identifier: str,
        name: str,
        positive: str,
        negative: str | None,
        *,
        level: str,
        parent: str | None = None,
        axis: str = "cell_type",
    ) -> None:
        labels.append(
            CorpusLabel.model_validate(
                {
                    "label_id": identifier,
                    "name": name,
                    "description": f"Synthetic definition for {name}; no biological claim.",
                    "source_id": "synthetic",
                    "source_label": name,
                    "species": ["Homo sapiens"], "tissues": ["blood"],
                    "positive": [positive],
                    "negative": [negative] if negative else [],
                    "level": level,
                    "parent_id": parent,
                    "axis": axis,
                    "exclusions": (exclusions or {}).get(identifier, []),
                    "status": "reviewed",
                    "reviewer": "synthetic-curator",
                    "rationale": "Known synthetic construction.",
                }
            )
        )

    for lineage, index in zip(["A", "B"], [0, 1], strict=True):
        label(
            f"lineage-{lineage}",
            f"Lineage {lineage}",
            f"lineage-{lineage}",
            f"lineage-{'B' if lineage == 'A' else 'A'}",
            level="L1",
        )
        for branch in range(2):
            branch_index = index * 2 + branch
            branch_id = f"branch-{branch_index}"
            label(
                branch_id,
                f"Branch {lineage}{branch}",
                branch_id,
                f"branch-{index * 2 + 1 - branch}",
                level="L2",
                parent=f"lineage-{lineage}",
            )
            for leaf in range(2):
                cluster = index * 4 + branch * 2 + leaf
                leaf_id = f"leaf-{cluster}"
                leaf_name = f"Leaf {lineage}{branch}{leaf}"
                label(
                    leaf_id,
                    leaf_name,
                    leaf_id,
                    f"leaf-{index * 4 + branch * 2 + 1 - leaf}",
                    level="L3",
                    parent=branch_id,
                )
                counts[
                    cluster * 20 : (cluster + 1) * 20, [index, 2 + branch_index, 6 + cluster]
                ] = [12, 9, 6]
                expected[str(cluster)] = [
                    f"Lineage {lineage}",
                    f"Branch {lineage}{branch}",
                    leaf_name,
                ]
    counts[:20, -1] = 4
    label("state-S", "S+", "state-S", None, level="L4", parent="leaf-0", axis="cell_state")
    source = tmp_path / "input.h5ad"
    cells = [f"cell-{i:03d}" for i in range(len(counts))]
    clusters = np.repeat([str(i) for i in range(8)], 20)
    if collapse_lineage:
        clusters[:80] = "A"
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame({"canonical_cluster": clusters, "lineage": "immune"}, index=cells),
        var=pd.DataFrame({"gene_symbol": features}, index=features),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    # The nearest raw cluster shares its L2 identity. A comparator must instead
    # test the nearest competing identity at the level under consideration.
    positions = np.repeat([0.0, 0.1, 10.0, 10.1, 100.0, 100.1, 110.0, 110.1], 20)
    adata.obsm["X_pca_harmony"] = np.column_stack([positions, np.zeros(len(counts))])
    adata.write_h5ad(source)
    corpus = CorpusSnapshot(
        corpus_id="synthetic-tree",
        revision="v1",
        sources=(
            CorpusSource(
                source_id="synthetic",
                name="Synthetic fixture",
                url="https://example.org/synthetic",
                version="1",
                sha256="a" * 64,
                retrieved_at="2026-11-07",
            ),
        ),
        labels=tuple(labels),
    )
    corpus_path = tmp_path / "corpus.json"
    publish_json(corpus_path, corpus.model_dump(mode="json"))
    config_path = _config(tmp_path, source, corpus_path)
    value = yaml.safe_load(config_path.read_text())
    value["knowledge"]["providers"] = [
        {
            "kind": "corpus",
            "path": str(corpus_path),
            "version": "v1",
            "sha256": sha256(corpus_path),
        }
    ]
    value["canonical_partitions"]["immune"]["expected_parent_count"] = (
        5 if collapse_lineage else 8
    )
    value["guidance"]["hierarchy_levels"] = ["L1", "L2", "L3"]
    value["guidance"]["context"]["tissue"] = "blood"  # decision row 16: corpus labels are tissue-scoped
    value["annotation"] = {
        "enabled": True,
        "assay": "transcriptome",
        "representation": "logcounts",
    }
    if not approved:
        value["guidance"].update({"interactive": True, "preauthorized": False})
    if classifier:
        model_path = tmp_path / "celltypist.pkl"
        model_path.write_bytes(b"synthetic model identity; predictions are mocked")
        value["annotation"]["classifiers"] = [
            {
                "name": "celltypist",
                "method": "celltypist",
                "model": str(model_path),
                "required": required,
                "label_map": {f"raw-{i}": expected[str(i)][-1] for i in range(8)},
            }
        ]

        def predictions(query, **_kwargs):
            if unavailable:
                return pd.DataFrame(), {
                    "status": "unavailable",
                    "reason": "synthetic missing dependency",
                }
            raw_labels = [
                "raw-7" if conflict else f"raw-{cluster}"
                for cluster in query.obs["canonical_cluster"]
            ]
            return pd.DataFrame(
                {
                    "barcode": query.obs_names.astype(str),
                    "predicted_label": raw_labels,
                    "cluster_id": query.obs["canonical_cluster"].astype(str).to_numpy(),
                }
            ), {"status": "complete"}

        monkeypatch.setattr("cell_curator.crosschecks.celltypist_hypothesis", predictions)
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    prepare_scene(config_path)
    context_path = tmp_path / "context-confirmation.json"
    context_path.write_text(json.dumps(value["guidance"]["context"]))
    confirm_context(config_path, input_path=context_path, reviewer="synthetic-reviewer")
    plan = construct_plan(config_path)
    assert plan["executable"] is True
    root = Path(value["run"]["output_root"]) / value["run"]["run_id"]
    return config_path, root, expected


def test_three_level_annotation_replays_and_keeps_parent_scopes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, expected = annotation_fixture(tmp_path, monkeypatch, classifier=True)
    source = Path(load_config(config_path).input.path)
    source_hash = sha256(source)
    run_classifier(config_path, "celltypist")
    for level in ("L1", "L2", "L3"):
        evidence = collect_level(config_path, level)
        resolve_level(config_path, level)
        assert all(item["status"] == "complete" for item in evidence["scopes"])
        scores = pd.read_csv(root / "annotation" / "levels" / level / "scores.tsv", sep="\t")
        for row in scores.itertuples():
            depth = int(level[1:]) - 1
            assert row.parent_scope == (
                "ROOT" if depth == 0 else expected[str(row.cluster_id)][depth - 1]
            )
        # Replaying already published stages must not rewrite evidence or decisions.
        digest = sha256(root / "annotation" / "levels" / level / "scores.tsv")
        assert collect_level(config_path, level) == evidence
        resolve_level(config_path, level)
        assert sha256(root / "annotation" / "levels" / level / "scores.tsv") == digest
    manifest = finalize_annotation(config_path)
    labels_path = root / "annotation" / "cell_labels.tsv"
    labels = pd.read_csv(labels_path, sep="\t", keep_default_na=False)
    assert len(labels) == 160
    for cluster, values in expected.items():
        block = labels[labels["cluster_id"].astype(str).eq(cluster)]
        assert block[["L1", "L2", "L3"]].drop_duplicates().values.tolist() == [values]
        assert set(block["L4"]) == {values[-1] + (" S+" if cluster == "0" else "")}
    assert set(labels["annotation_unit"]) == {"cluster_propagated_to_cells"}
    assert manifest["writeback_applied"] is False
    assert source_hash == sha256(source)
    assert finalize_annotation(config_path) == manifest


def test_corpus_and_classifier_resources_are_bound_before_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, _, _ = annotation_fixture(tmp_path, monkeypatch, classifier=True)
    model = load_config(config_path)
    contract = annotation_contract(model)
    assert contract["corpora"][0]["revision"] == "v1"
    Path(model.annotation.classifiers[0].model).write_bytes(b"changed model")
    with pytest.raises(ContractError, match="changed after planning"):
        run_classifier(config_path, "celltypist")


def test_gate_a_stops_annotation_before_collecting_label_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, _ = annotation_fixture(tmp_path, monkeypatch, approved=False)
    with pytest.raises(ContractError, match="approved|assumptions"):
        collect_level(config_path, "L1")
    assert not (root / "annotation" / "levels" / "L1" / "scores.tsv").exists()


def test_required_classifier_conflict_preserves_unknown(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, _ = annotation_fixture(
        tmp_path, monkeypatch, classifier=True, conflict=True
    )
    run_classifier(config_path, "celltypist")
    collect_level(config_path, "L1")
    resolve_level(config_path, "L1")
    labels = pd.read_csv(root / "annotation" / "levels" / "L1" / "labels.tsv", sep="\t")
    assert (
        labels[labels["cluster_id"].astype(str).isin(["0", "1", "2", "3"])]["label"]
        .eq("Unknown")
        .all()
    )
    scores = pd.read_csv(
        root / "annotation" / "levels" / "L1" / "scores.tsv", sep="\t", keep_default_na=False
    )
    supported = scores[
        scores["label"].eq("Lineage A") & scores["cluster_id"].astype(str).eq("0")
    ].iloc[0]
    assert "classifier:celltypist" in supported["gate_failures"]


def test_required_classifier_unavailable_stops_the_level(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, _ = annotation_fixture(
        tmp_path, monkeypatch, classifier=True, unavailable=True
    )
    with pytest.raises(ContractError, match="required classifier"):
        run_classifier(config_path, "celltypist")
    with pytest.raises(ContractError, match="required classifier"):
        collect_level(config_path, "L1")
    assert not (root / "annotation" / "levels" / "L1" / "scores.tsv").exists()


def test_required_classifier_must_run_before_marker_collection(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, _ = annotation_fixture(tmp_path, monkeypatch, classifier=True)
    with pytest.raises(ContractError, match="classifier.*missing"):
        collect_level(config_path, "L1")
    assert not (root / "annotation" / "levels" / "L1" / "scores.tsv").exists()


def test_missing_child_comparator_retains_the_supported_parent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, _ = annotation_fixture(tmp_path, monkeypatch, collapse_lineage=True)
    collect_level(config_path, "L1")
    resolve_level(config_path, "L1")
    evidence = collect_level(config_path, "L2")
    status = next(item for item in evidence["scopes"] if item["parent"] == "Lineage A")
    assert status["status"] == "unavailable_no_comparator"
    assert status["marker_methods"] == []
    resolve_level(config_path, "L2")
    collect_level(config_path, "L3")
    resolve_level(config_path, "L3")
    finalize_annotation(config_path)
    cells = pd.read_csv(root / "annotation" / "cell_labels.tsv", sep="\t")
    retained = cells[cells["cluster_id"].eq("A")]
    assert len(retained) == 80
    assert retained[["L1", "L2", "L3"]].eq("Lineage A").all().all()


def test_snakemake_executes_all_annotation_stages_and_resumes_without_work(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, expected = annotation_fixture(tmp_path, monkeypatch)
    work = tmp_path / "snakemake-work"
    work.mkdir()
    environment = {
        **os.environ,
        "PYTHONPATH": str(REPO_ROOT / "src"),
        "XDG_CACHE_HOME": str(tmp_path / "cache"),
        "PATH": f"{Path(sys.executable).parent}:{os.environ['PATH']}",
    }
    command = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(REPO_ROOT / "workflow" / "annotation.smk"),
        "--configfile",
        str(config_path),
        "--cores",
        "1",
        "--force-use-threads",
        "--directory",
        str(work),
    ]
    first = subprocess.run(
        command, env=environment, capture_output=True, text=True, timeout=360
    )
    assert first.returncode == 0, first.stdout + first.stderr
    cells = pd.read_csv(
        root / "annotation" / "cell_labels.tsv", sep="\t", keep_default_na=False
    )
    for cluster, values in expected.items():
        block = cells[cells["cluster_id"].astype(str).eq(cluster)]
        assert block[["L1", "L2", "L3"]].drop_duplicates().values.tolist() == [values]
        assert set(block["L4"]) == {values[-1] + (" S+" if cluster == "0" else "")}
    manifest = root / "annotation" / "annotation.manifest.json"
    digest = sha256(manifest)
    second = subprocess.run(
        command, env=environment, capture_output=True, text=True, timeout=120
    )
    assert second.returncode == 0, second.stdout + second.stderr
    assert "nothing to be done" in (second.stdout + second.stderr).lower()
    assert sha256(manifest) == digest


def test_advisory_classifier_conflict_is_a_flag_not_a_veto(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Decision row 10: a classifier that is not `required` cannot block a marker-supported call.
    config_path, root, _ = annotation_fixture(
        tmp_path, monkeypatch, classifier=True, conflict=True, required=False
    )
    run_classifier(config_path, "celltypist")
    collect_level(config_path, "L1")
    scores = pd.read_csv(
        root / "annotation" / "levels" / "L1" / "scores.tsv", sep="\t", keep_default_na=False
    )
    supported = scores[
        scores["label"].eq("Lineage A") & scores["cluster_id"].astype(str).eq("0")
    ].iloc[0]
    assert "classifier:celltypist" not in supported["gate_failures"]
    assert "classifier_disagreement:celltypist" in supported["gate_flags"]


def test_exclusion_is_tested_through_its_own_signature(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Decision row 30: a cluster that shows an excluded population's genes fails the label;
    # an exclusion whose genes the panel does not measure is flagged, not failed.
    config_path, root, _ = annotation_fixture(
        tmp_path,
        monkeypatch,
        exclusions={
            "lineage-A": [
                {"name": "self", "positive": ["lineage-A"]},
                {"name": "ghost", "positive": ["not-a-feature"]},
            ]
        },
    )
    collect_level(config_path, "L1")
    scores = pd.read_csv(
        root / "annotation" / "levels" / "L1" / "scores.tsv", sep="\t", keep_default_na=False
    )
    a0 = scores[scores["label"].eq("Lineage A") & scores["cluster_id"].astype(str).eq("0")].iloc[0]
    assert "excluded:self" in a0["gate_failures"]
    assert "excluded:ghost" not in a0["gate_failures"]
    assert "exclusion_unmeasured:ghost" in a0["gate_flags"]
    b4 = scores[scores["label"].eq("Lineage B") & scores["cluster_id"].astype(str).eq("4")].iloc[0]
    assert "excluded" not in b4["gate_failures"]
