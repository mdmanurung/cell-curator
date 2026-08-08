#!/usr/bin/env python3
"""Build or validate a lossless leaf mapping from individually accepted candidates."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from contracts import (
    OUTCOMES,
    ContractError,
    as_bool,
    atomic_write_json,
    atomic_write_tsv,
    read_config,
    require_columns,
    values_sha256,
)


def _split_ids(value: object) -> list[str]:
    return [token.strip() for token in str(value).split(";") if token.strip()]


def construct_mapping(
    *,
    canonical: pd.DataFrame,
    decisions: pd.DataFrame,
    registry: pd.DataFrame,
    membership: pd.DataFrame,
) -> pd.DataFrame:
    require_columns(canonical, {"barcode", "scope", "parent_id"}, artifact="canonical membership")
    require_columns(
        decisions,
        {"scope", "parent_id", "decision_outcome", "accepted_child_ids"},
        artifact="parent decisions",
    )
    require_columns(
        registry,
        {"scope", "candidate_key", "parent_id", "candidate_child_uid"},
        artifact="candidate registry",
    )
    require_columns(
        membership,
        {"scope", "candidate_key", "parent_id", "barcode", "candidate_membership"},
        artifact="candidate membership",
    )
    if canonical["barcode"].astype(str).duplicated().any():
        raise ContractError("canonical membership must contain one row per cell")
    if decisions.duplicated(["scope", "parent_id"]).any():
        raise ContractError("every parent must have exactly one decision")
    canonical_parents = set(zip(canonical["scope"].astype(str), canonical["parent_id"].astype(str)))
    decision_parents = set(zip(decisions["scope"].astype(str), decisions["parent_id"].astype(str)))
    if canonical_parents != decision_parents:
        raise ContractError("parent decisions do not cover canonical parents exactly")
    registry_by_uid = registry.set_index("candidate_child_uid", verify_integrity=True)
    rows: list[dict] = []
    for decision in decisions.itertuples(index=False):
        scope, parent = str(decision.scope), str(decision.parent_id)
        outcome = str(decision.decision_outcome)
        if outcome not in OUTCOMES:
            raise ContractError(f"invalid parent outcome: {outcome}")
        parent_frame = canonical[
            canonical["scope"].astype(str).eq(scope)
            & canonical["parent_id"].astype(str).eq(parent)
        ].copy()
        parent_uid = f"{scope}::{parent}"
        accepted = _split_ids(decision.accepted_child_ids)
        if outcome == "ACCEPT SPLIT" and not accepted:
            raise ContractError("ACCEPT SPLIT parent lacks accepted child IDs")
        if outcome != "ACCEPT SPLIT" and accepted:
            raise ContractError("unsplit parent declares accepted child IDs")
        assignments: dict[str, str] = {}
        if accepted:
            for child_uid in accepted:
                if child_uid not in registry_by_uid.index:
                    raise ContractError(f"accepted child is absent from registry: {child_uid}")
                candidate = registry_by_uid.loc[child_uid]
                if str(candidate["scope"]) != scope or str(candidate["parent_id"]) != parent:
                    raise ContractError("accepted child crosses its canonical parent")
                block = membership[
                    membership["candidate_key"].astype(str).eq(str(candidate["candidate_key"]))
                ]
                selected = block[block["candidate_membership"].map(as_bool)]
                for barcode in selected["barcode"].astype(str):
                    if barcode in assignments:
                        raise ContractError("accepted children overlap within a parent")
                    assignments[barcode] = child_uid
            remainder_uid = f"{parent_uid}::remainder"
        else:
            remainder_uid = parent_uid
        for barcode in parent_frame["barcode"].astype(str):
            leaf = assignments.get(barcode, remainder_uid)
            rows.append(
                {
                    "barcode": barcode,
                    "scope": scope,
                    "parent_id": parent,
                    "parent_uid": parent_uid,
                    "leaf_id": leaf,
                    "membership_role": (
                        "accepted_child"
                        if barcode in assignments
                        else ("parent_remainder" if accepted else "retained_parent")
                    ),
                    "decision_outcome": outcome,
                }
            )
    return pd.DataFrame(rows).sort_values("barcode", ignore_index=True)


def validate_mapping(canonical: pd.DataFrame, mapping: pd.DataFrame) -> dict:
    require_columns(
        mapping,
        {"barcode", "scope", "parent_id", "parent_uid", "leaf_id", "membership_role", "decision_outcome"},
        artifact="final mapping",
    )
    expected_barcodes = canonical["barcode"].astype(str)
    observed_barcodes = mapping["barcode"].astype(str)
    if observed_barcodes.duplicated().any() or set(observed_barcodes) != set(expected_barcodes):
        raise ContractError("mapping must contain exactly one row per input cell")
    expected = canonical.set_index(canonical["barcode"].astype(str))
    observed = mapping.set_index(mapping["barcode"].astype(str)).loc[expected.index]
    if not observed["scope"].astype(str).equals(expected["scope"].astype(str)) or not observed[
        "parent_id"
    ].astype(str).equals(expected["parent_id"].astype(str)):
        raise ContractError("mapping changed canonical parent membership")
    expected_parent_uid = observed["scope"].astype(str) + "::" + observed["parent_id"].astype(str)
    if not observed["parent_uid"].astype(str).equals(expected_parent_uid):
        raise ContractError("parent UID is inconsistent")
    canonical_counts = expected.groupby(["scope", "parent_id"], observed=True).size()
    mapped_counts = observed.groupby(["scope", "parent_id"], observed=True).size()
    if not canonical_counts.equals(mapped_counts):
        raise ContractError("parent counts are not conserved")
    for parent_uid, block in observed.groupby("parent_uid", observed=True):
        outcome = set(block["decision_outcome"].astype(str))
        if len(outcome) != 1:
            raise ContractError("a parent has multiple decision outcomes")
        if next(iter(outcome)) != "ACCEPT SPLIT":
            if set(block["leaf_id"].astype(str)) != {str(parent_uid)}:
                raise ContractError("unsplit parent membership changed")
        else:
            if not all(
                str(leaf).startswith(str(parent_uid) + "::")
                for leaf in block["leaf_id"].astype(str)
            ):
                raise ContractError("accepted split contains a cross-parent leaf")
    lineage_totals = (
        mapping.groupby("scope", observed=True)["leaf_id"].nunique().sort_index()
    )
    return {
        "schema_version": 1,
        "artifact_kind": "validated_final_leaf_mapping",
        "status": "complete",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_cells": len(mapping),
        "n_canonical_parents": int(mapping["parent_uid"].nunique()),
        "n_global_leaf_identities": int(mapping["leaf_id"].nunique()),
        "leaf_identities_by_scope": {
            str(key): int(value) for key, value in lineage_totals.items()
        },
        "ordered_barcode_sha256": values_sha256(list(mapping["barcode"].astype(str))),
        "sorted_barcode_sha256": values_sha256(
            list(mapping["barcode"].astype(str)), sort_values=True
        ),
        "invariants": {
            "one_row_per_input_cell": True,
            "canonical_parent_membership_unchanged": True,
            "parent_counts_conserved": True,
            "unsplit_membership_unchanged": True,
            "accepted_child_ids_drive_mapping": True,
            "failed_candidates_fold_to_remainder": True,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--candidate-membership", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--build-mapping", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        read_config(args.config)
        canonical = pd.read_csv(args.canonical, sep="\t", dtype=str, keep_default_na=False)
        decisions = pd.read_csv(args.decisions, sep="\t", dtype=str, keep_default_na=False)
        registry = pd.read_csv(args.registry, sep="\t", dtype=str, keep_default_na=False)
        membership = pd.read_csv(args.candidate_membership, sep="\t", dtype=str, keep_default_na=False)
        expected = construct_mapping(
            canonical=canonical,
            decisions=decisions,
            registry=registry,
            membership=membership,
        )
        if args.build_mapping:
            atomic_write_tsv(args.mapping, expected)
            mapping = expected
        else:
            mapping = pd.read_csv(args.mapping, sep="\t", dtype=str, keep_default_na=False)
            left = expected.sort_values("barcode", ignore_index=True).astype(str)
            right = mapping.sort_values("barcode", ignore_index=True)[left.columns].astype(str)
            if not left.equals(right):
                raise ContractError("provided mapping differs from accepted candidate decisions")
        atomic_write_json(args.output, validate_mapping(canonical, mapping))
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
