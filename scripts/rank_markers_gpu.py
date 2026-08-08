#!/usr/bin/env python3
"""Run native binary RAPIDS Wilcoxon and balanced logreg with no CPU fallback."""

from __future__ import annotations

import argparse
import gc
import gzip
import inspect
import os
import platform
import socket
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from contracts import ContractError, atomic_write_json, read_config, sha256, values_sha256


class GpuContractError(ContractError):
    """Raised before evidence creation when the GPU contract is invalid."""


def validate_runtime_contract_values(
    *,
    observed_version: str,
    expected_version: str,
    observed_distribution: str,
    expected_distribution: str,
    device_count: int,
    device_name: str,
    compute_capability: str,
    allowed_gpu_contracts: list[dict],
    observed_cuda_major: int,
    expected_cuda_major: int,
) -> dict:
    if not expected_version or observed_version != expected_version:
        raise GpuContractError(
            f"rapids-singlecell version mismatch: {observed_version!r} != {expected_version!r}"
        )
    if not expected_distribution or observed_distribution != expected_distribution:
        raise GpuContractError("RAPIDS distribution does not match configuration")
    if device_count != 1:
        raise GpuContractError(f"exactly one visible GPU is required, observed {device_count}")
    if observed_cuda_major != expected_cuda_major:
        raise GpuContractError("CUDA runtime major does not match configuration")
    if not allowed_gpu_contracts:
        raise GpuContractError("allowed_gpu_contracts must not be empty")
    for contract in allowed_gpu_contracts:
        name = str(contract["name_contains"])
        capability = str(contract["compute_capability"])
        if name.lower() in device_name.lower() and capability == compute_capability:
            return {
                "name_contains": name,
                "compute_capability": capability,
            }
    raise GpuContractError(
        f"GPU model/compute capability is not allowed: {device_name!r} {compute_capability}"
    )


def inspect_runtime(config: dict) -> tuple[object, object, dict]:
    marker = config["markers"]
    if marker.get("cpu_fallback") is not False:
        raise GpuContractError("markers.cpu_fallback must be false")
    import cupy as cp
    import rapids_singlecell as rsc

    observed_version = version("rapids-singlecell")
    distribution = str(marker["rapids_distribution"])
    observed_distribution = distribution
    device_count = int(cp.cuda.runtime.getDeviceCount())
    if device_count == 1:
        device = cp.cuda.Device(0)
        properties = cp.cuda.runtime.getDeviceProperties(0)
        raw_name = properties["name"]
        device_name = raw_name.decode() if isinstance(raw_name, bytes) else str(raw_name)
        capability = f"{device.compute_capability[0]}.{device.compute_capability[1]}"
    else:
        device_name, capability = "", ""
    runtime_version = int(cp.cuda.runtime.runtimeGetVersion())
    cuda_major = runtime_version // 1000
    validated = validate_runtime_contract_values(
        observed_version=observed_version,
        expected_version=str(marker["rapids_singlecell_version"]),
        observed_distribution=observed_distribution,
        expected_distribution=distribution,
        device_count=device_count,
        device_name=device_name,
        compute_capability=capability,
        allowed_gpu_contracts=list(marker["allowed_gpu_contracts"]),
        observed_cuda_major=cuda_major,
        expected_cuda_major=int(marker["cuda_major"]),
    )
    native = rsc.tl.rank_genes_groups
    return cp, rsc, {
        "device_count": device_count,
        "device_name": device_name,
        "compute_capability": capability,
        "validated_gpu_contract": validated,
        "rapids_singlecell": observed_version,
        "rapids_distribution": distribution,
        "cuda_major": cuda_major,
        "cuda_runtime_version": runtime_version,
        "cuda_driver_version": int(cp.cuda.runtime.driverGetVersion()),
        "native_api": f"{native.__module__}.{native.__name__}",
        "native_api_signature": str(inspect.signature(native)),
        "hostname": socket.gethostname(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
    }


def _result_frame(adata, key: str, *, method: str, metadata: dict) -> pd.DataFrame:
    result = adata.uns[key]
    names = result["names"]
    group = "0"
    if group not in (names.dtype.names or ()):
        raise GpuContractError("native result lacks candidate group 0")
    rows = []
    for rank, feature in enumerate(names[group], 1):
        feature = str(feature)
        def scalar(field: str) -> float:
            value = result.get(field)
            if value is None or group not in (value.dtype.names or ()):
                return float("nan")
            return float(value[group][rank - 1])
        rows.append(
            {
                **metadata,
                "method": method,
                "feature": feature,
                "rank": rank,
                "score": scalar("scores"),
                "logfoldchange": scalar("logfoldchanges") if method == "wilcoxon" else float("nan"),
                "pvalue": scalar("pvals"),
                "pvalue_adj": scalar("pvals_adj"),
                "effect_semantics": (
                    "native binary Wilcoxon candidate versus declared comparator"
                    if method == "wilcoxon"
                    else "balanced binary logistic-regression coefficient candidate versus declared comparator"
                ),
            }
        )
    return pd.DataFrame(rows)


def _write_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise GpuContractError(f"refusing to overwrite immutable output: {path}")
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    try:
        with gzip.open(partial, "wt", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, sep="\t", index=False)
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def run(args: argparse.Namespace) -> dict:
    outputs = (args.wilcoxon_output, args.logreg_output, args.receipt_output)
    if any(path.exists() for path in outputs):
        raise GpuContractError("refusing to overwrite GPU marker outputs")
    config = read_config(args.config)
    # Must be first runtime action: no input read and no output directory creation.
    cp, rsc, runtime = inspect_runtime(config)

    import anndata as ad

    adata = ad.read_h5ad(args.input)
    if args.group_key not in adata.obs:
        raise GpuContractError(f"input lacks binary group key {args.group_key!r}")
    groups = adata.obs[args.group_key].astype(str)
    if sorted(groups.unique()) != ["0", "1"]:
        raise GpuContractError("marker input must contain exactly binary groups 0 and 1")
    if args.matrix_source == "layer":
        matrix_key = args.matrix_key or args.representation
        if matrix_key not in adata.layers:
            raise GpuContractError("requested representation layer is absent")
        evidence_matrix = adata.layers[matrix_key]
    else:
        evidence_matrix = adata.X
    adata.obs[args.group_key] = pd.Categorical(groups, categories=["0", "1"], ordered=True)
    adata.X = evidence_matrix.copy()
    adata.layers.clear()
    adata.obsm.clear()
    adata.obsp.clear()
    gc.collect()
    marker = config["markers"]
    wilcoxon_key = "native_binary_wilcoxon"
    rsc.tl.rank_genes_groups(
        adata,
        groupby=args.group_key,
        groups=["0"],
        reference="1",
        n_genes=adata.n_vars,
        rankby_abs=False,
        pts=True,
        key_added=wilcoxon_key,
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        tie_correct=bool(marker["wilcoxon"]["tie_correct"]),
        use_continuity=bool(marker["wilcoxon"]["continuity_correct"]),
        use_raw=False,
        layer=None,
        chunk_size=int(marker["wilcoxon"]["chunk_size"]),
        multi_gpu=False,
    )
    logreg_key = "native_binary_logreg"
    rsc.tl.rank_genes_groups(
        adata,
        groupby=args.group_key,
        groups=["0"],
        reference="1",
        n_genes=adata.n_vars,
        rankby_abs=False,
        pts=True,
        key_added=logreg_key,
        method="logreg",
        use_raw=False,
        layer=None,
        class_weight="balanced",
        penalty="l2",
        C=float(marker["logreg"]["C"]),
        fit_intercept=True,
        max_iter=int(marker["logreg"]["retry_max_iter"]),
        tol=float(marker["logreg"]["tol"]),
        solver="qn",
    )
    metadata = {
        "lane_id": args.lane_id,
        "candidate_key": args.candidate_key,
        "comparator_role": args.comparator_role,
        "assay": args.assay,
        "assay_kind": args.assay_kind,
        "feature_space": args.feature_space,
        "representation": args.representation,
        "target_group": "0",
        "comparator_group": "1",
    }
    wilcoxon = _result_frame(adata, wilcoxon_key, method="wilcoxon", metadata=metadata)
    logreg = _result_frame(adata, logreg_key, method="logreg", metadata=metadata)
    _write_gzip(wilcoxon, args.wilcoxon_output)
    _write_gzip(logreg, args.logreg_output)
    receipt = {
        "schema_version": 1,
        "artifact_kind": "native_rapids_binary_marker_lane",
        "status": "complete",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **metadata,
        "methods": ["wilcoxon", "logreg"],
        "binary_groups": {"0": "candidate_child", "1": args.comparator_role},
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "ordered_barcode_sha256": values_sha256(list(adata.obs_names)),
        "input": {"path": str(args.input), "sha256": sha256(args.input)},
        "outputs": {
            "wilcoxon": {"path": str(args.wilcoxon_output), "sha256": sha256(args.wilcoxon_output)},
            "logreg": {"path": str(args.logreg_output), "sha256": sha256(args.logreg_output)},
        },
        "runtime": runtime,
        "software": {"python": platform.python_version()},
        "cpu_fallback_available": False,
        "cpu_fallback_used": False,
    }
    atomic_write_json(args.receipt_output, receipt)
    cp.get_default_memory_pool().free_all_blocks()
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--group-key", default="analysis_group")
    parser.add_argument("--lane-id", required=True)
    parser.add_argument("--candidate-key", required=True)
    parser.add_argument("--comparator-role", required=True, choices=("parent_remainder", "closest_sibling"))
    parser.add_argument("--assay", required=True)
    parser.add_argument("--assay-kind", required=True)
    parser.add_argument("--feature-space", required=True)
    parser.add_argument("--representation", required=True)
    parser.add_argument("--matrix-source", choices=("X", "layer", "x"), default="layer")
    parser.add_argument("--matrix-key", default="")
    parser.add_argument("--wilcoxon-output", required=True, type=Path)
    parser.add_argument("--logreg-output", required=True, type=Path)
    parser.add_argument("--receipt-output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except ContractError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
