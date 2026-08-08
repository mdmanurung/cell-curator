#!/usr/bin/env python3
"""Validate fail-closed RAPIDS receipts for configured RAPIDS evidence lanes."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from contracts import (
    COMPARATOR_ROLES,
    RAPIDS_BACKEND,
    ContractError,
    atomic_write_json,
    evidence_lanes,
    read_config,
    require_columns,
    sha256,
)


def validate_receipt(config: dict, lane: pd.Series, receipt: dict, root: Path) -> None:
    marker = config["markers"]
    expected = {
        "lane_id": str(lane["lane_id"]),
        "candidate_key": str(lane["candidate_key"]),
        "comparator_role": str(lane["comparator_role"]),
        "assay": str(lane["assay"]),
        "assay_kind": str(lane["assay_kind"]),
        "feature_space": str(lane["feature_space"]),
        "representation": str(lane["representation"]),
    }
    if (
        receipt.get("status") != "complete"
        or receipt.get("artifact_kind") != "native_rapids_binary_marker_lane"
    ):
        raise ContractError("RAPIDS receipt is incomplete or has wrong artifact kind")
    for key, value in expected.items():
        if str(receipt.get(key)) != value:
            raise ContractError(f"RAPIDS receipt field drifted: {key}")
    if receipt.get("methods") != ["wilcoxon", "logreg"]:
        raise ContractError("RAPIDS receipt lacks both native marker methods")
    if (
        receipt.get("cpu_fallback_available") is not False
        or receipt.get("cpu_fallback_used") is not False
    ):
        raise ContractError("RAPIDS receipt exposes a CPU fallback")
    runtime = receipt.get("runtime", {})
    if int(runtime.get("device_count", -1)) != 1:
        raise ContractError("RAPIDS receipt does not record exactly one device")
    if str(runtime.get("rapids_singlecell")) != str(
        marker["rapids_singlecell_version"]
    ):
        raise ContractError("RAPIDS version drifted")
    if str(runtime.get("rapids_distribution")) != str(
        marker["rapids_distribution"]
    ):
        raise ContractError("RAPIDS distribution drifted")
    if int(runtime.get("cuda_major", -1)) != int(marker["cuda_major"]):
        raise ContractError("CUDA major drifted")
    validated = runtime.get("validated_gpu_contract", {})
    allowed = {
        (str(item["name_contains"]), str(item["compute_capability"]))
        for item in marker["allowed_gpu_contracts"]
    }
    if (
        str(validated.get("name_contains")),
        str(validated.get("compute_capability")),
    ) not in allowed:
        raise ContractError("receipt device contract is not configured")
    for output in receipt.get("outputs", {}).values():
        path = Path(output["path"])
        if not path.is_absolute():
            path = root / path
        if not path.is_file() or sha256(path) != str(output["sha256"]):
            raise ContractError("RAPIDS marker output is missing or changed")


def validate_all(config: dict, lanes: pd.DataFrame, receipts_root: Path) -> dict:
    required = {
        "lane_id",
        "candidate_key",
        "comparator_role",
        "assay",
        "assay_kind",
        "feature_space",
        "representation",
        "evidence_lane_key",
        "backend",
        "receipt_path",
    }
    require_columns(lanes, required, artifact="lane registry")
    if lanes["lane_id"].duplicated().any():
        raise ContractError("lane registry duplicates lane IDs")

    configured = evidence_lanes(config)
    expected_lanes = {lane["lane_key"] for lane in configured}
    expected_roles = set(COMPARATOR_ROLES)
    for candidate, block in lanes.groupby("candidate_key", observed=True):
        observed = set(
            zip(
                block["comparator_role"].astype(str),
                block["evidence_lane_key"].astype(str),
            )
        )
        expected = {(role, lane) for role in expected_roles for lane in expected_lanes}
        if observed != expected or len(block) != len(expected):
            raise ContractError(
                f"candidate {candidate} lacks configured binary evidence lanes"
            )

    rapids = lanes[lanes["backend"].astype(str).eq(RAPIDS_BACKEND)]
    expected_rapids = {
        lane["lane_key"] for lane in configured if lane["backend"] == RAPIDS_BACKEND
    }
    for candidate, block in rapids.groupby("candidate_key", observed=True):
        observed = set(
            zip(
                block["comparator_role"].astype(str),
                block["evidence_lane_key"].astype(str),
            )
        )
        expected = {(role, lane) for role in expected_roles for lane in expected_rapids}
        if observed != expected:
            raise ContractError(f"candidate {candidate} lacks a required RAPIDS lane")
    for _, lane in rapids.iterrows():
        receipt_path = Path(str(lane["receipt_path"]))
        if not receipt_path.is_absolute():
            receipt_path = receipts_root / receipt_path
        if not receipt_path.is_file():
            raise ContractError(f"RAPIDS receipt is absent: {receipt_path}")
        validate_receipt(
            config, lane, json.loads(receipt_path.read_text()), receipts_root
        )
    return {
        "schema_version": 2,
        "artifact_kind": "validated_rapids_marker_receipts",
        "status": "complete",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_id": config["run"]["run_id"],
        "n_candidates": int(lanes["candidate_key"].nunique()),
        "n_configured_lanes": len(lanes),
        "n_rapids_lanes": len(rapids),
        "rapids_lanes_per_candidate": len(expected_rapids) * len(COMPARATOR_ROLES),
        "cpu_fallback_available": False,
        "cpu_fallback_used": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--lane-registry", required=True, type=Path)
    parser.add_argument("--receipts-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        config = read_config(args.config)
        lanes = pd.read_csv(
            args.lane_registry, sep="\t", dtype=str, keep_default_na=False
        )
        atomic_write_json(args.output, validate_all(config, lanes, args.receipts_root))
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
