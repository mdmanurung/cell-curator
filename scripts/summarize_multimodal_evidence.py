#!/usr/bin/env python3
"""Summarize absolute evidence for any configured assay representation."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from contracts import (
    ContractError,
    atomic_write_tsv,
    read_config,
    representation_decl,
)


def _dense_vector(matrix, column: int) -> np.ndarray:
    selected = matrix[:, column]
    if sparse.issparse(selected):
        return selected.toarray().reshape(-1).astype(float)
    return np.asarray(selected).reshape(-1).astype(float)


def _matrix(adata, declaration: dict):
    if str(declaration["source"]).lower() == "x":
        return adata.X
    key = str(declaration["key"])
    if key not in adata.layers:
        raise ContractError(f"declared layer is absent: {key}")
    return adata.layers[key]


def summarize_binary_features(
    *,
    expression,
    detection,
    groups: pd.Series,
    features: pd.Index,
    assay: str,
    assay_kind: str,
    feature_space: str,
    representation: str,
    detection_semantics: str = "unavailable",
) -> pd.DataFrame:
    labels = groups.astype(str).to_numpy()
    if sorted(set(labels)) != ["0", "1"]:
        raise ContractError("absolute evidence requires binary groups 0 and 1")
    if expression.shape[0] != len(labels):
        raise ContractError("expression and group dimensions differ")
    if detection is not None and expression.shape != detection.shape:
        raise ContractError("expression and detection dimensions differ")
    rows = []
    for column, feature in enumerate(features.astype(str)):
        values = _dense_vector(expression, column)
        detected_values = (
            None if detection is None else _dense_vector(detection, column)
        )
        if detected_values is not None and (detected_values < 0).any():
            raise ContractError("detection source must be nonnegative")
        summaries = {}
        for group in ("0", "1"):
            mask = labels == group
            summaries[group] = {
                "mean": float(values[mask].mean()),
                "median": float(np.median(values[mask])),
                "fraction_detected": (
                    float((detected_values[mask] > 0).mean())
                    if detected_values is not None
                    else float("nan")
                ),
            }
        rows.append(
            {
                "assay": assay,
                "assay_kind": assay_kind,
                "feature_space": feature_space,
                "representation": representation,
                "feature": feature,
                "target_mean": summaries["0"]["mean"],
                "comparator_mean": summaries["1"]["mean"],
                "delta_mean": summaries["0"]["mean"] - summaries["1"]["mean"],
                "target_median": summaries["0"]["median"],
                "comparator_median": summaries["1"]["median"],
                "target_fraction_detected": summaries["0"]["fraction_detected"],
                "comparator_fraction_detected": summaries["1"]["fraction_detected"],
                "delta_fraction_detected": summaries["0"]["fraction_detected"]
                - summaries["1"]["fraction_detected"],
                "detection_available": detection is not None,
                "detection_semantics": detection_semantics,
            }
        )
    return pd.DataFrame(rows)


def run(
    config: dict,
    input_path: Path,
    group_key: str,
    assay: str,
    representation: str,
) -> pd.DataFrame:
    import anndata as ad

    adata = ad.read_h5ad(input_path, backed="r")
    try:
        if group_key not in adata.obs:
            raise ContractError(f"input lacks group key {group_key!r}")
        declaration = representation_decl(config, assay, representation)
        expression = _matrix(adata, declaration)
        detection_decl = declaration.get("detection")
        detection = None
        semantics = "unavailable"
        if detection_decl is not None:
            detection = _matrix(adata, detection_decl)
            semantics = str(detection_decl.get("semantics", "> 0"))
        elif bool(declaration.get("detection_required", False)):
            raise ContractError("required detection source is absent")
        assay_decl = config["input"]["assays"][assay]
        return summarize_binary_features(
            expression=expression,
            detection=detection,
            groups=adata.obs[group_key],
            features=adata.var_names,
            assay=assay,
            assay_kind=str(assay_decl["kind"]),
            feature_space=str(
                declaration.get("feature_space", assay_decl["feature_space"])
            ),
            representation=representation,
            detection_semantics=semantics,
        )
    finally:
        adata.file.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--group-key", default="analysis_group")
    parser.add_argument("--assay", required=True)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        config = read_config(args.config)
        atomic_write_tsv(
            args.output,
            run(config, args.input, args.group_key, args.assay, args.representation),
        )
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
