#!/usr/bin/env python3
"""Build one label-free heterogeneity and technical audit row per parent."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .contracts import ContractError, atomic_write_tsv, read_config, require_columns


def normalized_entropy(values: pd.Series) -> float:
    counts = values.dropna().astype(str).value_counts()
    if counts.empty:
        return float("nan")
    if len(counts) == 1:
        return 0.0
    proportions = counts.to_numpy(dtype=float) / counts.sum()
    return -float(np.sum(proportions * np.log(proportions))) / math.log(len(counts))


def _top(values: pd.Series) -> tuple[str, float, float]:
    counts = values.dropna().astype(str).value_counts(normalize=True)
    if counts.empty:
        return "", float("nan"), float("nan")
    return (
        str(counts.index[0]),
        float(counts.iloc[0]),
        float(counts.iloc[1]) if len(counts) > 1 else 0.0,
    )


def _robust_z(values: pd.Series) -> pd.Series:
    median = float(values.median())
    mad = float((values - median).abs().median())
    if mad == 0:
        return pd.Series(0.0, index=values.index)
    return 0.67448975 * (values - median) / mad


def build_audit(config: dict, cells: pd.DataFrame) -> pd.DataFrame:
    keys = config["input"]["keys"]
    parent_key = keys["canonical_parent"]
    scope_key = keys["scope"]
    required = {
        "barcode",
        parent_key,
        scope_key,
    }
    optional_keys = {
        name: str(keys.get(name, "")).strip()
        for name in ("donor", "capture", "doublet_class", "sex", "batch")
    }
    required.update(value for value in optional_keys.values() if value)
    technical_metrics = {
        str(name): str(column)
        for name, column in config["input"].get("technical_metrics", {}).items()
        if str(column).strip()
    }
    required.update(technical_metrics.values())
    guides = list(config["input"].get("guide_columns", []))
    required.update(guides)
    require_columns(cells, required, artifact="cell metadata")
    if cells["barcode"].astype(str).duplicated().any():
        raise ContractError("cell metadata contains duplicate barcodes")

    thresholds = config["adaptive"]["audit"]
    rows: list[dict[str, Any]] = []
    for (scope, parent), frame in cells.groupby(
        [scope_key, parent_key], observed=True, sort=False
    ):
        donor_key = optional_keys["donor"]
        capture_key = optional_keys["capture"]
        doublet_key = optional_keys["doublet_class"]
        donors = (
            frame[donor_key].astype(str).value_counts(normalize=True)
            if donor_key
            else pd.Series(dtype=float)
        )
        captures = (
            frame[capture_key].astype(str).value_counts(normalize=True)
            if capture_key
            else pd.Series(dtype=float)
        )
        row: dict[str, Any] = {
            "scope": str(scope),
            "parent_id": str(parent),
            "parent_uid": f"{scope}::{parent}",
            "n_cells": len(frame),
            "n_donors": int(frame[donor_key].astype(str).nunique()) if donor_key else 0,
            "n_captures": int(frame[capture_key].astype(str).nunique()) if capture_key else 0,
            "max_donor_fraction": float(donors.iloc[0]) if not donors.empty else float("nan"),
            "max_capture_fraction": float(captures.iloc[0])
            if not captures.empty
            else float("nan"),
            "doublet_fraction": float(
                frame[doublet_key].astype(str).str.lower().eq("doublet").mean()
            )
            if doublet_key
            else float("nan"),
        }
        for metric, column in technical_metrics.items():
            row[f"median_technical__{metric}"] = float(frame[column].median())
        high_entropy = False
        meaningful_minority = False
        top_calls: set[str] = set()
        for guide in guides:
            top, top_fraction, second_fraction = _top(frame[guide])
            entropy = normalized_entropy(frame[guide])
            row[f"guide__{guide}__top"] = top
            row[f"guide__{guide}__top_fraction"] = top_fraction
            row[f"guide__{guide}__second_fraction"] = second_fraction
            row[f"guide__{guide}__entropy"] = entropy
            if top:
                top_calls.add(top)
            high_entropy |= bool(
                np.isfinite(entropy) and entropy >= float(thresholds["guide_entropy_min"])
            )
            meaningful_minority |= bool(
                np.isfinite(second_fraction)
                and second_fraction >= float(thresholds["meaningful_minority_fraction"])
            )
        row["guide_top_disagreement"] = len(top_calls) > 1
        row["high_guide_entropy"] = high_entropy
        row["meaningful_minority_identity"] = meaningful_minority
        row["doublet_enriched"] = bool(
            thresholds.get("doublet_fraction_max") is not None
            and np.isfinite(row["doublet_fraction"])
            and row["doublet_fraction"] >= float(thresholds["doublet_fraction_max"])
        )
        row["donor_dominated"] = bool(
            thresholds.get("donor_max_fraction") is not None
            and np.isfinite(row["max_donor_fraction"])
            and row["max_donor_fraction"] >= float(thresholds["donor_max_fraction"])
        )
        row["capture_dominated"] = bool(
            thresholds.get("capture_max_fraction") is not None
            and np.isfinite(row["max_capture_fraction"])
            and row["max_capture_fraction"] >= float(thresholds["capture_max_fraction"])
        )
        incompatible_columns = config["adaptive"]["audit"].get(
            "incompatible_program_columns", []
        )
        incompatible_fraction = max(
            [float(frame[column].astype(bool).mean()) for column in incompatible_columns]
            or [0.0]
        )
        row["incompatible_program_fraction"] = incompatible_fraction
        row["incompatible_program_enriched"] = incompatible_fraction >= float(
            config["adaptive"]["audit"].get("incompatible_program_fraction_min", 0.15)
        )
        sex_key = optional_keys["sex"]
        if sex_key and sex_key in frame:
            row["max_sex_fraction"] = float(
                frame[sex_key].astype(str).value_counts(normalize=True).iloc[0]
            )
        else:
            row["max_sex_fraction"] = float("nan")
        rows.append(row)

    result = pd.DataFrame(rows)
    if result.empty:
        raise ContractError("parent audit has no rows")
    median_metrics = [f"median_technical__{name}" for name in technical_metrics]
    for metric in median_metrics:
        result[f"{metric}_robust_z"] = _robust_z(result[metric])
    robust_metrics = [f"{metric}_robust_z" for metric in median_metrics]
    result["qc_outlier"] = (
        result[robust_metrics]
        .abs()
        .max(axis=1)
        .ge(float(thresholds.get("technical_metric_robust_z_max", 3.0)))
        if robust_metrics
        else False
    )
    reason_map = {
        "guide_top_disagreement": "guide_disagreement",
        "high_guide_entropy": "high_guide_entropy",
        "meaningful_minority_identity": "meaningful_minority_identity",
        "incompatible_program_enriched": "incompatible_programs",
        "doublet_enriched": "doublet_enrichment",
        "donor_dominated": "donor_dominance",
        "capture_dominated": "capture_dominance",
        "qc_outlier": "qc_or_complexity_outlier",
    }
    result["audit_flag_reasons"] = result.apply(
        lambda row: ";".join(
            reason for column, reason in reason_map.items() if bool(row[column])
        ),
        axis=1,
    )
    result["audit_flagged"] = result["audit_flag_reasons"].ne("")
    result["reviewed_outcome"] = ""
    result["outcome_status"] = np.where(
        result["audit_flagged"], "PENDING_EVIDENCE_REVIEW", "NOT_FLAGGED"
    )
    result.insert(0, "schema_version", 3)
    if result["parent_uid"].duplicated().any():
        raise ContractError("parent audit is not one row per canonical parent")
    return result.sort_values(["scope", "parent_id"], ignore_index=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        config = read_config(args.config)
        cells = pd.read_csv(args.cells, sep="\t")
        atomic_write_tsv(args.output, build_audit(config, cells))
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
