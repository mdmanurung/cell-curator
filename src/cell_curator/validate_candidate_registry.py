#!/usr/bin/env python3
"""Validate candidate routing, blank labels, membership, stability, and eligibility."""

from __future__ import annotations

import argparse
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .contracts import (
    ContractError,
    as_bool,
    atomic_write_json,
    read_config,
    require_blank_labels,
    require_columns,
    values_sha256,
)


def validate_registry(
    config: dict, audit: pd.DataFrame, registry: pd.DataFrame, membership: pd.DataFrame
) -> dict:
    require_columns(
        audit,
        {"scope", "parent_id", "audit_flagged"},
        artifact="parent audit",
    )
    require_columns(
        registry,
        {
            "scope",
            "candidate_key",
            "parent_id",
            "candidate_child_uid",
            "source_kind",
            "n_cells_in_parent",
            "n_parent_remainder",
            "n_donors",
            "n_captures",
            "stability_pass",
            "eligibility_counts_pass",
            "closest_sibling",
            "membership_sha256",
            "triage_action",
        },
        artifact="candidate registry",
    )
    require_columns(
        membership,
        {"scope", "candidate_key", "parent_id", "barcode", "candidate_membership"},
        artifact="candidate membership",
    )
    if registry["candidate_key"].duplicated().any():
        raise ContractError("candidate keys are not unique")
    require_blank_labels(registry)
    eligibility = config["adaptive"]["eligibility"]
    for row in registry.itertuples(index=False):
        if str(row.triage_action) != "FULL_TESTING":
            raise ContractError("registry may contain only FULL_TESTING candidates")
        if not as_bool(row.stability_pass) or not as_bool(row.eligibility_counts_pass):
            raise ContractError("registry candidate lacks stability or eligibility")
        if int(row.n_cells_in_parent) < int(eligibility["min_cells"]):
            raise ContractError("candidate cell count is below the configured minimum")
        if eligibility.get("min_donors") is not None and int(row.n_donors) < int(
            eligibility["min_donors"]
        ):
            raise ContractError("candidate donor count is below the configured minimum")
        if eligibility.get("min_captures") is not None and int(row.n_captures) < int(
            eligibility["min_captures"]
        ):
            raise ContractError("candidate capture count is below the configured minimum")
        if str(row.closest_sibling) == str(row.parent_id):
            raise ContractError("closest sibling equals the canonical parent")
        block = membership[membership["candidate_key"].astype(str).eq(str(row.candidate_key))]
        if block.empty or block["barcode"].astype(str).duplicated().any():
            raise ContractError("candidate membership is empty or duplicates parent cells")
        if set(block["scope"].astype(str)) != {str(row.scope)} or set(
            block["parent_id"].astype(str)
        ) != {str(row.parent_id)}:
            raise ContractError("candidate membership scope or parent drifted")
        selected = block[block["candidate_membership"].map(as_bool)]
        if len(selected) != int(row.n_cells_in_parent):
            raise ContractError("candidate membership count drifted from the registry")
        if len(block) - len(selected) != int(row.n_parent_remainder):
            raise ContractError("parent remainder count drifted from the registry")
        observed_hash = values_sha256(sorted(selected["barcode"].astype(str)))
        if observed_hash != str(row.membership_sha256):
            raise ContractError("candidate membership hash drifted")

    routed = set(
        zip(
            registry["scope"].astype(str),
            registry["parent_id"].astype(str),
            strict=True,
        )
    )
    flagged = {
        (str(row.scope), str(row.parent_id))
        for row in audit.itertuples(index=False)
        if as_bool(row.audit_flagged)
    }
    return {
        "schema_version": 3,
        "artifact_kind": "validated_candidate_registry",
        "status": "complete",
        "created_utc": datetime.now(UTC).isoformat(),
        "run_id": config["run"]["run_id"],
        "n_candidates": len(registry),
        "n_flagged_parents": len(flagged),
        "n_flagged_parents_with_candidates": len(flagged & routed),
        "invariants": {
            "candidate_keys_unique": True,
            "biological_labels_blank": True,
            "membership_one_row_per_parent_cell_per_candidate": True,
            "membership_hashes_valid": True,
            "stability_and_eligibility_pass": True,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        config = read_config(args.config)
        audit = pd.read_csv(args.audit, sep="\t", dtype=str, keep_default_na=False)
        registry = pd.read_csv(args.registry, sep="\t", dtype=str, keep_default_na=False)
        membership = pd.read_csv(args.membership, sep="\t", dtype=str, keep_default_na=False)
        atomic_write_json(args.output, validate_registry(config, audit, registry, membership))
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
