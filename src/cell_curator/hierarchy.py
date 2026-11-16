"""Strict parent-scoped hierarchy assignment, confidence, and durable label tables."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .contracts import ContractError
from .provenance import atomic_replace

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

