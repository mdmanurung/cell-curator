"""An object-oriented entry point over the configuration-driven workflow.

The functions in :mod:`cell_curator.api` each take a validated configuration
file, which is precise but makes the first five minutes with the package harder
than they need to be. :class:`CellCurator` builds that configuration from
keyword arguments and ordinary Python objects, then delegates every step to the
same functions the CLI and the Snakemake DAG call. It adds convenience, never
science: no method here computes evidence, relaxes a threshold, or decides a
label.

Two properties of the underlying design are deliberately preserved rather than
smoothed away:

* **The input is frozen.** An in-memory ``AnnData`` is written to the run
  directory and hashed before anything reads it, so a run always refers to one
  immutable object.
* **Labels stay gated.** There is no ``annotate_clusters()`` that returns
  finished labels. Evidence is computed for you; deciding what it means is the
  reviewer's job — a person, or the Claude/Codex agent driving the skill. That
  boundary is the reason the package exists, so the convenience layer does not
  cross it.

Nothing here talks to a third-party inference service. The annotator is whoever
reads the evidence.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .contracts import ContractError

TEMPLATE_RELATIVE = "share/cell-curator/assets/config.template.yaml"
SOURCE_RELATIVE = "skills/cell-curator/assets/config.template.yaml"


def default_config_template() -> dict[str, Any]:
    """Load the packaged strict configuration template.

    Resolution order is installed-distribution data first, then the source
    checkout, so an editable install and a wheel both work.
    """

    from importlib.metadata import PackageNotFoundError, distribution

    try:
        located = distribution("cell-curator").locate_file(TEMPLATE_RELATIVE)
    except PackageNotFoundError:
        located = None
    candidates = [Path(str(located))] if located is not None else []
    candidates.append(Path(__file__).resolve().parents[2] / SOURCE_RELATIVE)
    for candidate in candidates:
        if candidate.is_file():
            value = yaml.safe_load(candidate.read_text())
            if not isinstance(value, dict):
                raise ContractError(f"configuration template is not a mapping: {candidate}")
            return value
    raise ContractError(
        "could not locate assets/config.template.yaml in the installed "
        "distribution or the source checkout"
    )


def deep_update(base: dict[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively overlay ``overlay`` onto ``base`` and return ``base``."""

    for key, value in overlay.items():
        current = base.get(key)
        if isinstance(value, Mapping) and isinstance(current, dict):
            deep_update(current, value)
        else:
            base[key] = value
    return base


class CellCurator:
    """A configured, hash-frozen annotation run.

    Parameters mirror the questions the configuration actually asks. Anything
    not exposed here can still be set through ``overrides``, which is deep-merged
    onto the template last and therefore wins.

    Provide exactly one of ``adata`` or ``path``.
    """

    def __init__(
        self,
        *,
        adata: Any | None = None,
        path: str | Path | None = None,
        organism: str,
        tissue: str,
        run_id: str,
        cluster_key: str = "canonical_cluster",
        markers: str | Path | Mapping[str, Any] | None = None,
        output_root: str | Path = "results/cell_curator",
        mode: str = "evidence-only",
        scope_key: str | None = None,
        scope_name: str = "",
        embedding_key: str = "X_pca_harmony",
        layer: str = "logcounts",
        counts_layer: str = "counts",
        feature_id_column: str = "gene_symbol",
        developmental_stage: str = "",
        condition: str = "",
        technology: str = "dissociated",
        experimental_context: str = "",
        composition_covariates: Sequence[str] | None = None,
        reviewer: str = "",
        preauthorized: bool = False,
        assumptions: Sequence[str] = (),
        hierarchy_levels: Sequence[str] = ("L1",),
        writeback_prefix: str = "cell_curator",
        overrides: Mapping[str, Any] | None = None,
    ) -> None:
        if (adata is None) == (path is None):
            raise ContractError("provide exactly one of adata= or path=")
        self.run_id = run_id
        self.output_root = Path(output_root)
        self.workspace = self.output_root / run_id
        self.workspace.mkdir(parents=True, exist_ok=True)

        source, scope_key = self._materialize_input(
            adata=adata, path=path, scope_key=scope_key
        )
        marker_path = self._materialize_markers(markers)

        value = default_config_template()
        providers: list[dict[str, Any]] = []
        if marker_path is not None:
            providers.append(
                {
                    "kind": "local_markers",
                    "path": str(marker_path),
                    "version": f"{run_id}-markers",
                }
            )
        assay = value["input"]["assays"]["transcriptome"]
        representation = assay["representations"]["logcounts"]
        representation["key"] = layer
        representation["detection"]["key"] = counts_layer
        assay["feature_id_column"] = feature_id_column
        assay["graph_representations"] = [{"source": "obsm", "key": embedding_key}]

        deep_update(
            value,
            {
                "run": {
                    "run_id": run_id,
                    "output_root": str(self.output_root),
                    "mode": mode,
                },
                "knowledge": {"providers": providers},
                "input": {
                    "path": str(source),
                    "keys": {
                        "canonical_parent": cluster_key,
                        "scope": scope_key,
                        "embedding": embedding_key,
                    },
                },
                "canonical_partitions": {
                    scope_name or scope_key or "root": {
                        "parent_key": cluster_key,
                        "expected_parent_count": self._parent_count(source, cluster_key),
                    }
                },
                "guidance": {
                    "interactive": not preauthorized,
                    "preauthorized": preauthorized,
                    "reviewer": reviewer,
                    "pause_after_each_level": not preauthorized,
                    "assumptions": list(assumptions),
                    "hierarchy_levels": list(hierarchy_levels),
                    "context": {
                        "organism": organism,
                        "tissue": tissue,
                        "developmental_stage": developmental_stage,
                        "condition": condition,
                        "technology": technology,
                        "experimental_context": experimental_context,
                        "composition_covariates": (
                            list(composition_covariates)
                            if composition_covariates is not None
                            else [scope_key]
                        ),
                        "vocabulary_contract": (
                            "local signed marker programs" if marker_path else ""
                        ),
                        "reference_contract": "local_only",
                        "visualization_key": embedding_key,
                        "writeback_prefix": writeback_prefix,
                    },
                },
            },
        )
        value["canonical_partitions"].pop("example_scope", None)
        if overrides:
            deep_update(value, overrides)

        self.config_path = self.workspace / "config.yaml"
        self.config_path.write_text(yaml.safe_dump(value, sort_keys=False))
        self.validate()

    # -- construction helpers -------------------------------------------------

    DEFAULT_SCOPE_KEY = "cell_curator_scope"

    def _materialize_input(
        self, *, adata: Any | None, path: str | Path | None, scope_key: str | None
    ) -> tuple[Path, str]:
        """Write the input to the run directory so the run can hash-freeze it.

        ``input.keys.scope`` names the obs column that partitions the frozen
        parents. Most datasets do not carry one, so when ``scope_key`` is not
        given a single-scope column is added *to the run directory's copy*. The
        caller's object and file are never modified.
        """

        destination = self.workspace / f"{self.run_id}.input.h5ad"
        if path is not None:
            resolved = Path(path)
            if not resolved.is_file():
                raise ContractError(f"input does not exist: {resolved}")
            if scope_key:
                self._require_obs_column(resolved, scope_key)
                return resolved.resolve(), scope_key
            import anndata as ad

            adata = ad.read_h5ad(resolved)
        assert adata is not None  # guaranteed by the exactly-one check in __init__

        if destination.exists():
            raise ContractError(
                f"{destination} already exists; a run identifier is single-use because "
                "the input is hash-frozen. Choose a new run_id, or pass path= with an "
                "explicit scope_key= to read an existing object in place."
            )
        if scope_key:
            if scope_key not in adata.obs:
                raise ContractError(
                    f"scope_key {scope_key!r} is absent from obs; "
                    f"available keys include {sorted(adata.obs.columns)[:12]}"
                )
            adata.write_h5ad(destination)
            return destination.resolve(), scope_key

        resolved_scope = self.DEFAULT_SCOPE_KEY
        if resolved_scope in adata.obs:
            raise ContractError(
                f"obs already carries {resolved_scope!r}; pass scope_key={resolved_scope!r} "
                "to use it, or scope_key= naming a different partition column"
            )
        # Categorical because composition covariates must be categorical.
        adata.obs[resolved_scope] = pd.Categorical(["all"] * adata.n_obs)
        try:
            adata.write_h5ad(destination)
        finally:
            # Restore the caller's object: only the run directory copy carries it.
            del adata.obs[resolved_scope]
        return destination.resolve(), resolved_scope

    @staticmethod
    def _require_obs_column(source: Path, column: str) -> None:
        import anndata as ad

        observations = ad.read_h5ad(source, backed="r").obs
        if column not in observations:
            raise ContractError(
                f"scope_key {column!r} is absent from obs; "
                f"available keys include {sorted(observations.columns)[:12]}"
            )

    def _materialize_markers(
        self, markers: str | Path | Mapping[str, Any] | None
    ) -> Path | None:
        if markers is None:
            return None
        if isinstance(markers, Mapping):
            destination = self.workspace / f"{self.run_id}.markers.json"
            destination.write_text(json.dumps(markers, indent=2, sort_keys=True))
            return destination.resolve()
        resolved = Path(markers)
        if not resolved.is_file():
            raise ContractError(f"marker programs do not exist: {resolved}")
        return resolved.resolve()

    @staticmethod
    def _parent_count(source: Path, cluster_key: str) -> int:
        """Count observed clusters so the declared partition matches the data."""

        import anndata as ad

        observations = ad.read_h5ad(source, backed="r").obs
        if cluster_key not in observations:
            raise ContractError(
                f"cluster_key {cluster_key!r} is absent from obs; "
                f"available keys include {sorted(observations.columns)[:12]}"
            )
        return int(observations[cluster_key].astype(str).nunique())

    # -- workflow -------------------------------------------------------------

    def validate(self) -> Any:
        """Load and strictly validate the generated configuration."""

        from .config import load_config

        return load_config(self.config_path)

    @property
    def run_directory(self) -> Path:
        """The immutable run root that every artifact is published under."""

        from .pipeline import run_root

        return run_root(self.validate())

    def inspect(self) -> dict[str, Any]:
        """Report structure, keys, and representations without changing anything."""

        from .api import inspect_input

        return inspect_input(self.config_path)

    def freeze(self) -> dict[str, Any]:
        """Create the hash-frozen working copy, membership artifacts, and QC."""

        from .api import canonicalize_input, compute_input_qc

        canonicalized = canonicalize_input(self.config_path)
        canonicalized["qc"] = compute_input_qc(self.config_path)
        return canonicalized

    def scene(self) -> dict[str, Any]:
        """Persist Stage 1 context, compute detection, and the assumptions ledger."""

        from .api import prepare_annotation_scene

        return prepare_annotation_scene(self.config_path)

    def plan(self) -> dict[str, Any]:
        """Persist Stage 2 vocabulary, evidence routes, and open assumptions."""

        from .api import plan_annotation

        return plan_annotation(self.config_path)

    def prepare(self) -> dict[str, Any]:
        """Freeze the input, then run the scene and plan stages in order.

        Evidence generation requires both stages, so this removes the one piece
        of ordering a caller would otherwise have to know. An interactive run
        (the default) still stops here until the context is confirmed and the
        assumptions gate is signed.
        """

        frozen = self.freeze()
        scene = self.scene()
        plan = self.plan()
        return {"frozen": frozen, "scene": scene, "plan": plan}

    def approve_assumptions(
        self, *, reviewer: str, scope: str = "*", preauthorized: bool = False
    ) -> dict[str, Any]:
        """Sign the assumptions gate. Required before any biological label."""

        from .api import approve_annotation_assumptions

        return approve_annotation_assumptions(
            self.config_path,
            reviewer=reviewer,
            scope=scope,
            preauthorized=preauthorized,
        )

    def markers(self, **kwargs: Any) -> Any:
        """Run the signed marker profile and configured cross-check surface."""

        from .api import run_marker_tools

        return run_marker_tools(self.config_path, **kwargs)

    # -- adaptive refinement --------------------------------------------------

    def audit_parents(self) -> Any:
        """Flag clusters that look impure, before proposing any split.

        Returns one row per frozen parent cluster carrying the signals that
        justify looking closer: donor and capture dominance, doublet fraction,
        incompatible-program fraction, technical-metric outliers, and guide
        entropy. A flag is a reason to investigate, never a decision — the
        thresholds live under ``adaptive.audit`` in the configuration.
        """

        from .pipeline import audit

        return audit(self.config_path)

    def propose_subclusters(self) -> tuple[Any, Any, Any]:
        """Discover bounded subcluster candidates for flagged parents.

        Canonical parents stay frozen. This first tries to reuse a stored
        higher-resolution clustering whose children match a parent cleanly
        (purity, fraction, and adjacency criteria under
        ``adaptive.stored_candidates``); only when no stored family matches does
        it fall back to a bounded parent-local reclustering pass
        (``adaptive.parent_local``, one pass by default, rejecting one-cluster
        outcomes and requiring stability above the configured threshold).

        Returns ``(registry, candidates, evaluations)`` and writes the
        comparator manifest each candidate is later tested against. Nothing is
        accepted here — every candidate must still clear
        ``evidence.required_candidate_gates`` and reviewer sign-off.
        """

        from .pipeline import refine

        return refine(self.config_path)

    def candidate_evidence(self, lane_id: str) -> Path:
        """Run one candidate-versus-comparator evidence lane by manifest id."""

        from .pipeline import candidate_evidence_lane

        return candidate_evidence_lane(self.config_path, comparison_lane_id=lane_id)

    def decide_subclusters(self) -> tuple[Any, Any]:
        """Derive split/retain outcomes for every candidate from frozen evidence.

        Each parent resolves to one of ``ACCEPT SPLIT``, ``RETAIN PARENT``,
        ``DOUBLET/MIXED``, or ``TECHNICAL/UNRESOLVED``. A parent is split only
        when the evidence gates pass; impurity that turns out to be technical or
        doublet-driven is reported as such rather than carved into subtypes.
        """

        from .api import propose_annotations

        return propose_annotations(self.config_path)

    def refine_clusters(self) -> dict[str, Any]:
        """Audit, propose, score, and decide subclusters in one call.

        Convenience over :meth:`audit_parents`, :meth:`propose_subclusters`,
        :meth:`evidence`, and :meth:`decide_subclusters`. It reports outcomes;
        it does not write labels or modify the input.
        """

        audited = self.audit_parents()
        registry, candidates, evaluations = self.propose_subclusters()
        self.evidence()
        decisions, parents = self.decide_subclusters()
        return {
            "audited_parents": audited,
            "candidate_registry": registry,
            "candidates": candidates,
            "evaluations": evaluations,
            "decisions": decisions,
            "parent_outcomes": parents,
        }

    def evidence(self) -> tuple[Any, Any]:
        """Run every executable evidence lane and return the consolidated tables."""

        from .api import generate_evidence

        return generate_evidence(self.config_path)

    def review_packet(self, *, strict: bool = True) -> Path:
        """Build evidence cards, critic packets, and the self-contained report."""

        from .api import build_review_artifacts

        return build_review_artifacts(self.config_path, strict=strict)

    def status(self) -> dict[str, Any]:
        """Return the current stage, blockers, and the next action to take."""

        from .api import annotation_progress

        return annotation_progress(self.config_path)

    def validate_run(self) -> dict[str, Any]:
        """Check completeness, hashes, mapping conservation, and review state."""

        from .api import validate_run

        return validate_run(self.config_path)

    def preview_writeback(self, *args: Any, **kwargs: Any) -> Any:
        """Return the source-preserving write-back diff without applying it."""

        from .api import preview_writeback

        return preview_writeback(self.config_path, *args, **kwargs)

    def write_back(self, *args: Any, **kwargs: Any) -> Any:
        """Apply an explicitly approved write-back to a new object by default."""

        from .api import apply_writeback

        return apply_writeback(self.config_path, *args, **kwargs)

    def __repr__(self) -> str:
        return f"CellCurator(run_id={self.run_id!r}, config={str(self.config_path)!r})"
