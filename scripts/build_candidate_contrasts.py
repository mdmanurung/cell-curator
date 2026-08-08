#!/usr/bin/env python3
"""Build immutable binary candidate-vs-remainder and candidate-vs-sibling lanes."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from contracts import (
    COMPARATOR_ROLES,
    ContractError,
    as_bool,
    atomic_write_tsv,
    evidence_lanes,
    read_config,
    require_columns,
    values_sha256,
)


def build_binary_groups(
    *,
    cells: pd.DataFrame,
    membership: pd.DataFrame,
    candidate: pd.Series,
    comparator_role: str,
) -> pd.DataFrame:
    if comparator_role not in COMPARATOR_ROLES:
        raise ContractError(f"unsupported comparator role: {comparator_role}")
    scope = str(candidate["scope"])
    parent = str(candidate["parent_id"])
    block = membership[
        membership["candidate_key"].astype(str).eq(str(candidate["candidate_key"]))
    ].copy()
    if block.empty or block["barcode"].astype(str).duplicated().any():
        raise ContractError("candidate membership is empty or duplicated")
    is_candidate = block["candidate_membership"].map(as_bool)
    target = set(block.loc[is_candidate, "barcode"].astype(str))
    remainder = set(block.loc[~is_candidate, "barcode"].astype(str))
    if not target or not remainder:
        raise ContractError("candidate or parent remainder is empty")
    scope_cells = cells[cells["scope"].astype(str).eq(scope)]
    if comparator_role == "parent_remainder":
        comparator = remainder
    else:
        sibling = str(candidate["closest_sibling"])
        if not sibling or sibling == parent:
            raise ContractError("closest sibling is absent or equals parent")
        comparator = set(
            scope_cells.loc[
                scope_cells["parent_id"].astype(str).eq(sibling), "barcode"
            ].astype(str)
        )
    if not comparator or target & comparator:
        raise ContractError("binary comparator is empty or overlaps the candidate")
    selected = scope_cells[
        scope_cells["barcode"].astype(str).isin(target | comparator)
    ].copy()
    selected["analysis_group"] = selected["barcode"].astype(str).map(
        lambda barcode: "0" if barcode in target else "1"
    )
    selected["analysis_role"] = selected["analysis_group"].map(
        {"0": "candidate_child", "1": comparator_role}
    )
    if set(selected["analysis_group"]) != {"0", "1"}:
        raise ContractError("binary comparison lacks group 0 or group 1")
    if selected["barcode"].astype(str).duplicated().any():
        raise ContractError("binary comparison duplicates barcodes")
    return selected.sort_values("barcode", ignore_index=True)


def build_contrasts(
    config: dict,
    cells: pd.DataFrame,
    registry: pd.DataFrame,
    membership: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    require_columns(cells, {"barcode", "scope", "parent_id"}, artifact="cells")
    require_columns(
        registry,
        {"scope", "candidate_key", "parent_id", "closest_sibling"},
        artifact="candidate registry",
    )
    require_columns(
        membership,
        {"candidate_key", "barcode", "candidate_membership"},
        artifact="candidate membership",
    )
    if cells["barcode"].astype(str).duplicated().any():
        raise ContractError("cells must have exactly one row per input barcode")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ContractError(f"refusing to write into nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    lanes = evidence_lanes(config)
    rows: list[dict] = []
    for _, candidate in registry.iterrows():
        for comparator_role in COMPARATOR_ROLES:
            groups = build_binary_groups(
                cells=cells,
                membership=membership,
                candidate=candidate,
                comparator_role=comparator_role,
            )
            base = output_dir / str(candidate["candidate_key"]) / comparator_role
            groups_path = base / "groups.tsv"
            atomic_write_tsv(groups_path, groups)
            target_hash = values_sha256(
                sorted(
                    groups.loc[
                        groups["analysis_group"].eq("0"), "barcode"
                    ].astype(str)
                )
            )
            comparator_hash = values_sha256(
                sorted(
                    groups.loc[
                        groups["analysis_group"].eq("1"), "barcode"
                    ].astype(str)
                )
            )
            for lane in lanes:
                assay = str(lane["assay"])
                rep = str(lane["representation"])
                lane_id = (
                    f"{candidate['candidate_key']}__{comparator_role}__{lane['lane_key']}"
                )
                rows.append(
                    {
                        "lane_id": lane_id,
                        "scope": str(candidate["scope"]),
                        "candidate_key": str(candidate["candidate_key"]),
                        "parent_id": str(candidate["parent_id"]),
                        "comparator_role": comparator_role,
                        "assay": assay,
                        "assay_kind": str(lane["assay_kind"]),
                        "feature_space": str(lane["feature_space"]),
                        "representation": rep,
                        "evidence_lane_key": str(lane["lane_key"]),
                        "backend": str(lane["backend"]),
                        "required_for_acceptance": bool(
                            lane["required_for_acceptance"]
                        ),
                        "ontology_compatible": bool(lane["ontology_compatible"]),
                        "detection_required": bool(lane["detection_required"]),
                        "groups_path": str(groups_path),
                        "target_n_cells": int(groups["analysis_group"].eq("0").sum()),
                        "comparator_n_cells": int(groups["analysis_group"].eq("1").sum()),
                        "target_barcode_sha256": target_hash,
                        "comparator_barcode_sha256": comparator_hash,
                        "receipt_path": f"{lane_id}/receipt.json",
                        "wilcoxon_path": f"{lane_id}/wilcoxon.tsv.gz",
                        "logreg_path": f"{lane_id}/logreg.tsv.gz",
                    }
                )
    result = pd.DataFrame(rows)
    expected = len(registry) * len(COMPARATOR_ROLES) * len(lanes)
    if len(result) != expected or result["lane_id"].duplicated().any():
        raise ContractError("binary comparator lane registry is incomplete or duplicated")
    atomic_write_tsv(output_dir / "comparator_manifest.tsv", result)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        config = read_config(args.config)
        cells = pd.read_csv(args.cells, sep="\t", dtype=str, keep_default_na=False)
        registry = pd.read_csv(
            args.registry, sep="\t", dtype=str, keep_default_na=False
        )
        membership = pd.read_csv(
            args.membership, sep="\t", dtype=str, keep_default_na=False
        )
        build_contrasts(config, cells, registry, membership, args.output_dir)
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
