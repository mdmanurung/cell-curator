"""Strict parent-scoped hierarchy assignment, confidence, and durable label tables."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .contracts import ContractError, identifier_text
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


def validate_parent_vocabulary(vocabulary: pd.DataFrame) -> None:
    required = {"level", "parent_scope", "label"}
    missing = required.difference(vocabulary.columns)
    if missing:
        raise ContractError(f"hierarchy vocabulary misses columns: {sorted(missing)}")
    values = vocabulary.loc[:, ["level", "parent_scope", "label"]].astype(str)
    if values.apply(lambda column: column.str.strip().eq("").any()).any():
        raise ContractError("hierarchy vocabulary contains blank level/parent/label values")
    if "axis" in vocabulary.columns:
        # This vocabulary assigns identity, so a cell state program reaching it would
        # become a cell's type -- the one direction the separation exists to prevent.
        # A frame without the column carries no state rows for it to admit.
        state = sorted(
            set(
                vocabulary.loc[~vocabulary["axis"].astype(str).eq("cell_type"), "label"].astype(
                    str
                )
            )
        )
        if state:
            raise ContractError(
                "hierarchy vocabulary assigns cell identity, so every row must carry "
                f"axis cell_type; cell state programs present: {state}"
            )
    parents_per_label = values.groupby("label", observed=True)["parent_scope"].nunique()
    invalid = sorted(parents_per_label[parents_per_label > 1].index.astype(str))
    if invalid:
        raise ContractError(f"a child label may belong to only one parent scope: {invalid}")


def identity_vocabulary(vocabulary: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Restrict a planned vocabulary to identity rows, refusing a cross-axis parent.

    The planned vocabulary records every marker program the providers declared,
    cell state included, because the plan is that record. Assignment needs only the
    identity rows, and dropping the rest is safe unless a cell type program descends
    from a cell state one: filtering would then leave the survivor without its parent.
    Nothing constrains ``MarkerProgram.parent`` to one axis, so that link is reachable
    and is refused here rather than silently producing a broken hierarchy.
    """

    if "axis" not in vocabulary.columns:
        return vocabulary, []
    is_identity = vocabulary["axis"].astype(str).eq("cell_type")
    excluded = sorted(set(vocabulary.loc[~is_identity, "label"].astype(str)))
    if not excluded:
        return vocabulary, []
    identity = vocabulary.loc[is_identity].reset_index(drop=True)
    dropped = set(excluded)
    orphaned = sorted(
        f"{label} under {parent}"
        for label, parent in zip(
            identity["label"].astype(str),
            identity["parent_scope"].astype(str),
            strict=True,
        )
        if parent in dropped
    )
    if orphaned:
        raise ContractError(
            "a cell type program descends from a cell state program, so identity and "
            f"state cannot be separated: {orphaned}"
        )
    return identity, excluded


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
    # Two parent scopes are not identities: an abstained parent ("Unknown"), whose
    # accepted children resolve what it could not name, and the label itself, when a
    # split keeps its parent's identity. Neither gives a label a second parent.
    scope = non_unknown["parent_scope"].astype(str)
    lineage = non_unknown[~scope.eq("Unknown") & ~scope.eq(non_unknown["label"].astype(str))]
    parents = lineage.groupby("label", observed=True)["parent_scope"].nunique()
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
                "cl_id": identifier_text(proposal.get("cl_id")),
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
            from .data import resolve_obsm

            spatial, _, _ = resolve_obsm(loaded, spatial_key)
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
    "labels_from_existing_proposals",
    "profile_refined_clusters",
    "scaffold_labels_table",
    "validate_hierarchy_nesting",
    "identity_vocabulary",
    "validate_parent_vocabulary",
]
