"""Strict parent-scoped hierarchy assignment, confidence, and durable label tables."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .contracts import ContractError
from .provenance import atomic_publish, atomic_replace, publish_json

LABEL_COLUMNS = (
    "level",
    "parent_scope",
    "cluster_id",
    "label",
    "cell_state",
    "cl_id",
    "confidence",
    "margin",
    "specificity",
    "agreement",
    "coverage",
    "rationale",
    "rejected_alternatives",
    "reviewer",
)


def advisory_confidence(
    *,
    margin: float,
    specificity: float,
    agreement: float,
    coverage: float,
) -> float:
    """Combine four evidence dimensions without pretending to be calibrated probability."""

    margin_component = 1.0 - np.exp(-max(float(margin), 0.0))
    values = np.clip(
        [margin_component, specificity, agreement, coverage],
        0.0,
        1.0,
    )
    return float(np.dot(values, [0.35, 0.25, 0.25, 0.15]))


def classify_parallel_state(
    *,
    mixed_score: float = 0.0,
    technical_effect: float = 0.0,
    depth_robust_z: float = 0.0,
    mitochondrial_robust_z: float = 0.0,
    cycling_fraction: float = 0.0,
    stress_score: float = 0.0,
    gradient_score: float = 0.0,
) -> dict[str, Any]:
    """Classify technical/state axes independently from biological identity."""

    flags: list[str] = []
    if mixed_score >= 0.6:
        flags.append("mixed/doublet")
    if technical_effect >= 0.3:
        flags.append("technical")
    if depth_robust_z <= -3.0 or mitochondrial_robust_z >= 3.0:
        flags.append("low-quality")
    if cycling_fraction >= 0.5:
        flags.append("cycling")
    if stress_score >= 0.6:
        flags.append("stressed")
    if gradient_score >= 0.6:
        flags.append("gradient")
    state = flags[0] if flags else "normal"
    identity_must_abstain = state in {"mixed/doublet", "technical", "low-quality"}
    return {
        "cell_state": state,
        "state_flags": flags,
        "identity_must_abstain": identity_must_abstain,
    }


def profile_ambiguity(
    cell_scores: pd.DataFrame,
    *,
    first_label: str,
    second_label: str,
) -> dict[str, Any]:
    """Distinguish mixed co-support from a continuous/uncertain score gradient."""

    required = {first_label, second_label}
    missing = required.difference(cell_scores.columns)
    if missing:
        raise ContractError(f"ambiguity profile misses score columns: {sorted(missing)}")
    first = cell_scores[first_label].astype(float).to_numpy()
    second = cell_scores[second_label].astype(float).to_numpy()
    if not len(first):
        raise ContractError("ambiguity profile requires at least one cell")
    high_first = first >= np.nanquantile(first, 0.75)
    high_second = second >= np.nanquantile(second, 0.75)
    co_high = float(np.mean(high_first & high_second))
    contrast = first - second
    middle = float(np.mean(np.abs(contrast) <= max(float(np.nanstd(contrast)) * 0.25, 1e-8)))
    correlation = float(np.corrcoef(first, second)[0, 1]) if len(first) > 1 else 0.0
    correlation = 0.0 if not np.isfinite(correlation) else correlation
    if co_high >= 0.2:
        interpretation = "mixed_coexpression"
    elif middle >= 0.35 or abs(correlation) >= 0.7:
        interpretation = "continuous_gradient"
    else:
        interpretation = "separable_or_sparse"
    return {
        "interpretation": interpretation,
        "co_high_fraction": co_high,
        "middle_contrast_fraction": middle,
        "score_correlation": correlation,
        "profiled_before_unknown": True,
    }


def validate_parent_vocabulary(vocabulary: pd.DataFrame) -> None:
    required = {"level", "parent_scope", "label"}
    missing = required.difference(vocabulary.columns)
    if missing:
        raise ContractError(f"hierarchy vocabulary misses columns: {sorted(missing)}")
    values = vocabulary.loc[:, ["level", "parent_scope", "label"]].astype(str)
    if values.apply(lambda column: column.str.strip().eq("").any()).any():
        raise ContractError("hierarchy vocabulary contains blank level/parent/label values")
    parents_per_label = values.groupby("label", observed=True)["parent_scope"].nunique()
    invalid = sorted(parents_per_label[parents_per_label > 1].index.astype(str))
    if invalid:
        raise ContractError(f"a child label may belong to only one parent scope: {invalid}")


def assign_parent_scoped_level(
    scores: pd.DataFrame,
    vocabulary: pd.DataFrame,
    *,
    level: str,
    parent_scope: str,
    minimum_confidence: float = 0.55,
    minimum_margin: float = 0.1,
    state_by_cluster: dict[str, dict[str, Any]] | None = None,
) -> pd.DataFrame:
    """Assign one level while strictly filtering candidates to the declared parent."""

    validate_parent_vocabulary(vocabulary)
    required = {"cluster_id", "label", "score"}
    missing = required.difference(scores.columns)
    if missing:
        raise ContractError(f"hierarchy score table misses columns: {sorted(missing)}")
    if "level" in scores.columns:
        observed_levels = set(scores["level"].dropna().astype(str))
        if observed_levels and observed_levels != {level}:
            raise ContractError(
                "hierarchy score table declares a different level than its execution "
                f"scope: expected {level!r}, observed {sorted(observed_levels)}"
            )
    if "parent_scope" in scores.columns:
        observed_parents = set(scores["parent_scope"].dropna().astype(str))
        if observed_parents and observed_parents != {parent_scope}:
            raise ContractError(
                "hierarchy score table declares a different parent than its execution "
                f"scope: expected {parent_scope!r}, observed {sorted(observed_parents)}"
            )
    candidates = vocabulary[
        vocabulary["level"].astype(str).eq(level)
        & vocabulary["parent_scope"].astype(str).eq(parent_scope)
    ]
    if candidates.empty:
        raise ContractError(
            f"hierarchy vocabulary has no candidates for {level} within {parent_scope}"
        )
    permitted = set(candidates["label"].astype(str))
    block = scores[scores["label"].astype(str).isin(permitted)].copy()
    if block.empty:
        raise ContractError(
            f"no executable evidence remains after parent scoping for {level}::{parent_scope}"
        )
    rows: list[dict[str, Any]] = []
    state_by_cluster = state_by_cluster or {}
    cl_lookup = dict(
        zip(
            candidates["label"].astype(str),
            candidates.get("cl_id", pd.Series("", index=candidates.index)).astype(str),
            strict=True,
        )
    )
    for cluster_id, cluster_scores in block.groupby("cluster_id", observed=True, sort=True):
        ordered = cluster_scores.sort_values(["score", "label"], ascending=[False, True])
        winner = ordered.iloc[0]
        runner = ordered.iloc[1] if len(ordered) > 1 else None
        margin = (
            float(winner["score"] - runner["score"])
            if runner is not None
            else float(winner["score"])
        )
        specificity = float(winner.get("specificity", winner.get("de_overlap", 0.0)))
        agreement = float(winner.get("agreement", winner.get("crosscheck_agreement", 0.0)))
        coverage = float(winner.get("coverage", winner.get("positive_coverage", 0.0)))
        confidence = advisory_confidence(
            margin=margin,
            specificity=specificity,
            agreement=agreement,
            coverage=coverage,
        )
        state = state_by_cluster.get(str(cluster_id), {"cell_state": "normal"})
        rejected = [
            f"{row.label} (score={float(row.score):.3f})"
            for row in ordered.iloc[1:4].itertuples(index=False)
        ]
        abstain_reasons = []
        if confidence < minimum_confidence:
            abstain_reasons.append("advisory confidence below threshold")
        if margin < minimum_margin:
            abstain_reasons.append("evidence margin below threshold")
        if bool(state.get("identity_must_abstain", False)):
            abstain_reasons.append(f"terminal parallel state: {state['cell_state']}")
        label = "Unknown" if abstain_reasons else str(winner["label"])
        rationale = (
            "; ".join(abstain_reasons)
            if abstain_reasons
            else (
                f"Parent-scoped winner {label}; margin={margin:.3f}, specificity={specificity:.3f}, "
                f"agreement={agreement:.3f}, coverage={coverage:.3f}."
            )
        )
        rows.append(
            {
                "level": level,
                "parent_scope": parent_scope,
                "cluster_id": str(cluster_id),
                "label": label,
                "cell_state": str(state.get("cell_state", "normal")),
                "cl_id": "" if label == "Unknown" else cl_lookup.get(label, ""),
                "confidence": confidence,
                "margin": margin,
                "specificity": specificity,
                "agreement": agreement,
                "coverage": coverage,
                "rationale": rationale,
                "rejected_alternatives": " | ".join(rejected) or "none measured",
                "reviewer": "",
            }
        )
    return pd.DataFrame(rows, columns=LABEL_COLUMNS)


def scaffold_labels_table(
    path: str | Path,
    proposals: pd.DataFrame,
    *,
    retain_unmentioned: bool = True,
) -> pd.DataFrame:
    """Update the living labels source of truth without losing sibling scopes.

    Existing rows outside the supplied scope set are retained by default so that
    independently executed parents can accumulate in one durable level table.  A
    named reviewer marks a matching row as manually curated; its biological call
    and reasoning are then preserved across automatic rebuilds.
    """

    destination = Path(path)
    missing = set(LABEL_COLUMNS).difference(proposals.columns)
    if missing:
        raise ContractError(f"label proposals miss durable columns: {sorted(missing)}")
    frame = proposals.loc[:, list(LABEL_COLUMNS)].copy()
    keys = ["level", "parent_scope", "cluster_id"]
    for column in keys:
        frame[column] = frame[column].astype(str)
        if frame[column].str.strip().eq("").any():
            raise ContractError(f"label proposals contain blank durable key {column!r}")
    if frame.duplicated(keys).any():
        raise ContractError("label proposals contain duplicate level/parent/cluster rows")
    if destination.is_file():
        existing = pd.read_csv(destination, sep="\t", keep_default_na=False)
        if set(LABEL_COLUMNS).difference(existing.columns):
            raise ContractError("existing labels table does not satisfy the durable schema")
        existing = existing.loc[:, list(LABEL_COLUMNS)].copy()
        for column in keys:
            existing[column] = existing[column].astype(str)
            if existing[column].str.strip().eq("").any():
                raise ContractError(
                    f"existing labels table contains blank durable key {column!r}"
                )
        if existing.duplicated(keys).any():
            duplicate = existing.loc[existing.duplicated(keys, keep=False), keys].iloc[0]
            key = tuple(str(duplicate[column]) for column in keys)
            raise ContractError(f"existing labels table duplicates key {key}")
        indexed = existing.set_index(keys)
        manually_curated = (
            "label",
            "cell_state",
            "cl_id",
            "confidence",
            "rationale",
            "rejected_alternatives",
            "reviewer",
        )
        for row_index, row in frame.iterrows():
            key = tuple(str(row[item]) for item in keys)
            if key not in indexed.index:
                continue
            old = indexed.loc[key]
            if isinstance(old, pd.DataFrame):
                raise ContractError(f"existing labels table duplicates key {key}")
            reviewer = str(old["reviewer"]).strip()
            if reviewer and reviewer.upper() not in {"TODO", "PENDING"}:
                for column in manually_curated:
                    value = str(old[column]).strip()
                    if value and value.upper() not in {"TODO", "PENDING"}:
                        frame.at[row_index, column] = old[column]
        if retain_unmentioned:
            incoming_keys = pd.MultiIndex.from_frame(frame[keys])
            existing_keys = pd.MultiIndex.from_frame(existing[keys])
            untouched = existing.loc[~existing_keys.isin(incoming_keys)]
            frame = pd.concat([untouched, frame], ignore_index=True)
    frame = frame.sort_values(keys, kind="stable").reset_index(drop=True)
    atomic_replace(destination, frame.to_csv(sep="\t", index=False).encode())
    return frame


def recursive_parent_scoped_assignment(
    *,
    vocabulary: pd.DataFrame,
    score_tables: dict[tuple[str, str], pd.DataFrame],
    state_tables: dict[tuple[str, str], dict[str, dict[str, Any]]] | None = None,
    minimum_confidence: float = 0.55,
    minimum_margin: float = 0.1,
) -> dict[tuple[str, str], pd.DataFrame]:
    """Execute an arbitrary-depth hierarchy from independently checkpointable scopes."""

    validate_parent_vocabulary(vocabulary)
    outputs: dict[tuple[str, str], pd.DataFrame] = {}
    state_tables = state_tables or {}
    pending: list[tuple[str, str]] = [("L1", "ROOT")]
    visited: set[tuple[str, str]] = set()
    while pending:
        scope = pending.pop(0)
        if scope in visited:
            continue
        visited.add(scope)
        if scope not in score_tables:
            # A parent without an executable score table is a successful retained leaf,
            # not permission to score it against some other parent's candidates.
            continue
        level, parent = scope
        assigned = assign_parent_scoped_level(
            score_tables[scope],
            vocabulary,
            level=level,
            parent_scope=parent,
            minimum_confidence=minimum_confidence,
            minimum_margin=minimum_margin,
            state_by_cluster=state_tables.get(scope),
        )
        outputs[scope] = assigned
        try:
            level_number = int(level.removeprefix("L"))
        except ValueError as exc:
            raise ContractError(f"hierarchy level must use L<number>: {level}") from exc
        next_level = f"L{level_number + 1}"
        available_parents = set(
            vocabulary.loc[
                vocabulary["level"].astype(str).eq(next_level), "parent_scope"
            ].astype(str)
        )
        for label in sorted(set(assigned["label"]).intersection(available_parents)):
            pending.append((next_level, label))
    return outputs


def validate_hierarchy_nesting(
    level_assignments: Iterable[pd.DataFrame],
    *,
    cell_membership: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Enforce one parent per label and optional one-cell-one-leaf membership."""

    frames = [frame.copy() for frame in level_assignments]
    if not frames:
        raise ContractError("hierarchy validation requires at least one labels table")
    combined = pd.concat(frames, ignore_index=True)
    required = {"level", "parent_scope", "label", "cluster_id"}
    if required.difference(combined.columns):
        raise ContractError("hierarchy labels lack level/parent/label/cluster columns")
    non_unknown = combined[~combined["label"].astype(str).eq("Unknown")]
    parents = non_unknown.groupby("label", observed=True)["parent_scope"].nunique()
    invalid = sorted(parents[parents > 1].index.astype(str))
    if invalid:
        raise ContractError(f"strict hierarchy nesting failed for labels: {invalid}")
    invariants = {
        "one_parent_per_label": True,
        "unknown_permitted": True,
        "one_cell_one_leaf": None,
    }
    if cell_membership is not None:
        required_cells = {"barcode", "leaf_id"}
        if required_cells.difference(cell_membership.columns):
            raise ContractError("cell membership requires barcode and leaf_id")
        if cell_membership["barcode"].astype(str).duplicated().any():
            raise ContractError("one-cell-one-leaf validation found duplicate barcodes")
        if cell_membership["leaf_id"].astype(str).str.strip().eq("").any():
            raise ContractError("one-cell-one-leaf validation found blank leaves")
        invariants["one_cell_one_leaf"] = True
    return {
        "schema_version": 3,
        "artifact_kind": "cell_curator_hierarchy_validation",
        "status": "complete",
        "n_levels": int(combined["level"].nunique()),
        "n_labels": int(non_unknown["label"].nunique()),
        "invariants": invariants,
    }


def labels_from_existing_proposals(
    proposals: pd.DataFrame,
    *,
    level: str = "L1",
    parent_scope: str = "ROOT",
) -> pd.DataFrame:
    """Bridge the existing multimodal proposal abstraction into durable labels."""

    rows = []
    for proposal in proposals.to_dict("records"):
        alternatives = proposal.get("rejected_alternatives", "")
        if isinstance(alternatives, (list, tuple)):
            alternatives = " | ".join(map(str, alternatives))
        rows.append(
            {
                "level": level,
                "parent_scope": parent_scope,
                "cluster_id": str(proposal["cluster_id"]),
                "label": str(proposal.get("cell_type", "Unknown")),
                "cell_state": str(proposal.get("cell_state", "normal")),
                "cl_id": str(proposal.get("cl_id", "")),
                "confidence": float(proposal.get("confidence", 0.0)),
                "margin": float(proposal.get("margin", 0.0)),
                "specificity": float(proposal.get("specificity", 0.0)),
                "agreement": float(
                    proposal.get(
                        "agreement", 1.0 - float(proposal.get("normalized_entropy", 1.0))
                    )
                ),
                "coverage": float(proposal.get("positive_coverage", 0.0)),
                "rationale": str(proposal.get("rationale", "")),
                "rejected_alternatives": str(alternatives or "none measured"),
                "reviewer": str(proposal.get("reviewer", "")),
            }
        )
    return pd.DataFrame(rows, columns=LABEL_COLUMNS)


def profile_refined_clusters(config_path: str | Path) -> dict[str, Any]:
    """Build durable per-child QC, composition, and spatial summaries."""

    from .config import load_config
    from .contracts import sha256
    from .data import close_input, load_input

    config = load_config(config_path)
    root = Path(config.run.output_root) / config.run.run_id
    registry_path = root / "discovery" / "candidate_registry.tsv"
    membership_path = root / "discovery" / "candidate_membership.tsv"
    if not registry_path.is_file() or not membership_path.is_file():
        from .pipeline import refine

        refine(Path(config_path))
    registry = pd.read_csv(registry_path, sep="\t", keep_default_na=False)
    membership = pd.read_csv(membership_path, sep="\t", keep_default_na=False)
    qc_path = root / "qc" / "cell_metrics.tsv"
    qc = (
        pd.read_csv(qc_path, sep="\t", keep_default_na=False).set_index("barcode")
        if qc_path.is_file()
        else pd.DataFrame()
    )
    context_path = root / "context.json"
    covariates = list(config.guidance.context.composition_covariates)
    if context_path.is_file():
        scene = json.loads(context_path.read_text())
        covariates = list(
            scene.get("biological_context", {}).get("composition_covariates", covariates)
        )
    loaded = load_input(config, backed=True)
    try:
        obs = loaded.root.obs.copy()
        obs.index = obs.index.astype(str)
        spatial_key = config.input.keys.spatial
        spatial: np.ndarray | None = None
        if spatial_key:
            if spatial_key in loaded.root.obsm:
                spatial = np.asarray(loaded.root.obsm[spatial_key])
            else:
                spatial = next(
                    (
                        np.asarray(adata.obsm[spatial_key])
                        for adata in loaded.assays.values()
                        if spatial_key in adata.obsm
                    ),
                    None,
                )
        spatial_frame = (
            pd.DataFrame(
                spatial,
                index=obs.index,
                columns=[f"spatial_{index + 1}" for index in range(spatial.shape[1])],
            )
            if spatial is not None
            else pd.DataFrame(index=obs.index)
        )
        rows: list[dict[str, Any]] = []
        for candidate in registry.to_dict("records"):
            candidate_key = str(candidate["candidate_key"])
            block = membership[
                membership["candidate_key"].astype(str).eq(candidate_key)
                & membership["candidate_membership"]
                .astype(str)
                .str.casefold()
                .isin({"true", "1"})
            ]
            barcodes = block["barcode"].astype(str).tolist()
            identity = {
                "scope": str(candidate["scope"]),
                "parent_id": str(candidate["parent_id"]),
                "candidate_key": candidate_key,
                "candidate_child_uid": str(candidate["candidate_child_uid"]),
                "n_cells": len(barcodes),
            }
            rows.append(
                {
                    **identity,
                    "panel": "size",
                    "metric": "n_cells",
                    "category": "",
                    "value": float(len(barcodes)),
                }
            )
            if not qc.empty:
                available = [barcode for barcode in barcodes if barcode in qc.index]
                numeric = qc.loc[available].select_dtypes(include=[np.number])
                for column, value in numeric.median(axis=0).items():
                    rows.append(
                        {
                            **identity,
                            "panel": "qc",
                            "metric": str(column),
                            "category": "median",
                            "value": float(value),
                        }
                    )
            for covariate in covariates:
                if covariate not in obs:
                    continue
                available = [barcode for barcode in barcodes if barcode in obs.index]
                fractions = (
                    obs.loc[available, covariate].astype(str).value_counts(normalize=True)
                )
                for category, value in fractions.items():
                    rows.append(
                        {
                            **identity,
                            "panel": "composition",
                            "metric": covariate,
                            "category": str(category),
                            "value": float(value),
                        }
                    )
            if not spatial_frame.empty:
                available = [barcode for barcode in barcodes if barcode in spatial_frame.index]
                for dimension, value in spatial_frame.loc[available].mean(axis=0).items():
                    rows.append(
                        {
                            **identity,
                            "panel": "spatial",
                            "metric": str(dimension),
                            "category": "centroid",
                            "value": float(value),
                        }
                    )
    finally:
        close_input(loaded)
    columns = [
        "scope",
        "parent_id",
        "candidate_key",
        "candidate_child_uid",
        "n_cells",
        "panel",
        "metric",
        "category",
        "value",
    ]
    profiles = pd.DataFrame(rows, columns=columns)
    output = root / "discovery" / "cluster_profiles.tsv"
    atomic_publish(output, profiles.to_csv(sep="\t", index=False).encode())
    marker_receipts = sorted((root / "evidence" / "receipts" / "candidates").glob("*.json"))
    manifest = {
        "schema_version": 3,
        "artifact_kind": "cell_curator_refined_cluster_profiles",
        "status": "complete",
        "n_candidates": int(registry["candidate_key"].nunique()) if not registry.empty else 0,
        "n_rows": len(profiles),
        "panels": sorted(profiles["panel"].unique()) if not profiles.empty else [],
        "profile_path": str(output),
        "profile_sha256": sha256(output),
        "candidate_marker_receipts": [str(path) for path in marker_receipts],
        "limitations": (
            []
            if marker_receipts
            else [
                "per-child marker receipts are not available until candidate evidence lanes run"
            ]
        ),
    }
    publish_json(root / "discovery" / "cluster_profiles.manifest.json", manifest)
    return manifest


__all__ = [
    "LABEL_COLUMNS",
    "advisory_confidence",
    "assign_parent_scoped_level",
    "classify_parallel_state",
    "labels_from_existing_proposals",
    "profile_ambiguity",
    "profile_refined_clusters",
    "recursive_parent_scoped_assignment",
    "scaffold_labels_table",
    "validate_hierarchy_nesting",
    "validate_parent_vocabulary",
]
