"""Deterministic classifier diagnostics alongside marker-derived labels.

Predictions remain hypotheses. Diagnostics do not change marker labels, turn
classifier vote shares into probabilities, or assume independent methods.
"""

from __future__ import annotations

import json
from collections.abc import Mapping

import numpy as np
import pandas as pd

from .contracts import ContractError


def _ancestors(label: str, parents: Mapping[str, str]) -> list[str]:
    result = [label]
    while result[-1] in parents:
        parent = parents[result[-1]]
        if not parent or parent in result:
            raise ContractError("classifier reconciliation parent links contain a cycle or blank parent")
        result.append(parent)
    return result


def reconcile_predictions(
    marker_labels: pd.DataFrame,
    predictions: Mapping[str, pd.DataFrame],
    *,
    label_mapping: Mapping[str, Mapping[str, str]] | None = None,
    parent_by_label: Mapping[str, str] | None = None,
    cell_membership: pd.DataFrame | None = None,
    min_support: float = 0.6,
) -> pd.DataFrame:
    """Attach mapped classifier composition and conflicts without overriding calls.

    ``label_mapping[provider][raw_label]`` must contain reviewed exact mappings.
    Unmapped predictions stay unresolved. When supplied, cell_membership binds
    every classifier to exactly the declared barcodes and cluster membership.
    A shared ancestor is reported as a suggestion, never silently assigned.
    """
    if not np.isfinite(min_support) or not 0 < min_support <= 1:
        raise ContractError("classifier min_support must be finite and in (0, 1]")
    if {"cluster_id", "label"} - set(marker_labels):
        raise ContractError("marker labels require cluster_id and label")
    if marker_labels.empty or marker_labels[["cluster_id", "label"]].isna().any().any():
        raise ContractError("marker labels must be nonempty with explicit cluster and label values")
    labels = marker_labels.copy()
    labels["cluster_id"] = labels["cluster_id"].astype(str)
    if labels["cluster_id"].duplicated().any():
        raise ContractError("classifier reconciliation requires one marker label per cluster")
    mappings = label_mapping or {}
    parents = parent_by_label or {}
    for child in parents:
        _ancestors(child, parents)
    known = set(labels["label"].astype(str)) | set(parents) | set(parents.values())
    known.discard("Unknown")
    for mapping in mappings.values():
        if any(not str(raw).strip() or target not in known for raw, target in mapping.items()):
            raise ContractError("classifier mapping contains a blank alias or unknown canonical label")
    membership: pd.Series | None = None
    if cell_membership is not None:
        if {"barcode", "cluster_id"} - set(cell_membership):
            raise ContractError("cell membership requires barcode and cluster_id")
        if cell_membership["barcode"].astype(str).duplicated().any():
            raise ContractError("cell membership barcodes must be unique")
        membership = cell_membership.set_index(cell_membership["barcode"].astype(str))["cluster_id"].astype(str)
        if set(membership) != set(labels["cluster_id"]):
            raise ContractError("cell membership clusters differ from marker label clusters")
    normalized: dict[str, pd.DataFrame] = {}
    for provider, frame in sorted(predictions.items()):
        if {"barcode", "predicted_label"} - set(frame):
            raise ContractError(f"{provider} predictions require barcode and predicted_label")
        block = frame.copy()
        block["barcode"] = block["barcode"].astype(str)
        if block["barcode"].duplicated().any() or block["barcode"].str.strip().eq("").any():
            raise ContractError(f"{provider} prediction barcodes must be nonblank and unique")
        if membership is not None:
            if set(block["barcode"]) != set(membership.index):
                raise ContractError(f"{provider} prediction barcodes differ from cell membership")
            expected = block["barcode"].map(membership)
            if "cluster_id" in block and not block["cluster_id"].astype(str).eq(expected).all():
                raise ContractError(f"{provider} prediction cluster membership differs")
            block["cluster_id"] = expected
        if "cluster_id" not in block or set(block["cluster_id"].astype(str)) != set(labels["cluster_id"]):
            raise ContractError(f"{provider} prediction clusters differ from marker label clusters")
        if block["predicted_label"].isna().any() or block["predicted_label"].astype(str).str.strip().eq("").any():
            raise ContractError(f"{provider} predictions contain missing labels; use explicit Unknown")
        block["cluster_id"] = block["cluster_id"].astype(str)
        raw = block["predicted_label"].astype(str)
        mapping = mappings.get(provider, {})
        block["canonical_label"] = raw.map(mapping).fillna(raw.where(raw.isin(known), "Unknown"))
        normalized[provider] = block
    rows: list[dict[str, object]] = []
    for row in labels.itertuples(index=False):
        cluster = str(row.cluster_id)
        marker_label = str(row.label)
        summaries = {}
        confident = []
        conflicts = []
        support = []
        dependent = []
        for provider, block in normalized.items():
            selected = block.loc[block["cluster_id"].eq(cluster)]
            fractions = selected["canonical_label"].value_counts(normalize=True)
            ordered = sorted(((str(k), float(v)) for k, v in fractions.items()), key=lambda item: (-item[1], item[0]))
            winner, fraction = ordered[0]
            tied = len(ordered) > 1 and ordered[1][1] == fraction
            summaries[provider] = {"composition": dict(ordered), "n_cells": len(selected), "unresolved": winner == "Unknown" or tied or fraction < min_support}
            if "popv_celltypist_prediction" in selected:
                dependent.append("popv includes CELLTYPIST; no independence from celltypist assumed")
            if summaries[provider]["unresolved"]:
                continue
            confident.append(winner)
            compatible = marker_label in _ancestors(winner, parents) or winner in _ancestors(marker_label, parents)
            if marker_label != "Unknown" and compatible:
                support.append(provider)
            elif marker_label != "Unknown":
                conflicts.append(provider)
        common = ""
        if confident:
            paths = [_ancestors(label, parents) for label in confident]
            common = next((label for label in paths[0] if all(label in path for path in paths[1:])), "")
        rows.append({
            "classifier_summary": json.dumps(summaries, sort_keys=True),
            "classifier_supporting_providers": "|".join(support),
            "classifier_conflicting_providers": "|".join(conflicts),
            "classifier_conflict": bool(conflicts),
            "classifier_parent_suggestion": common,
            "classifier_dependency_notes": "|".join(dependent),
            "classifier_evidence_role": "corroboration_not_independent_replication",
        })
    return pd.concat([labels.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
