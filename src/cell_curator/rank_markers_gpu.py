#!/usr/bin/env python3
"""Validate the native RAPIDS runtime contract with no CPU fallback."""

from __future__ import annotations

import inspect
import os
import socket
from importlib.metadata import distribution, version
from typing import Any

from .contracts import ContractError


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


def inspect_runtime(config: dict[str, Any]) -> tuple[Any, Any, dict[str, Any]]:
    marker = config["markers"]
    if marker.get("cpu_fallback") is not False:
        raise GpuContractError("markers.cpu_fallback must be false")
    import cupy as cp
    import rapids_singlecell as rsc

    observed_version = version("rapids-singlecell")
    expected_distribution = str(marker["rapids_distribution"])
    package = distribution("rapids-singlecell")
    observed_distribution = str(package.metadata["Name"])
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
        expected_distribution=expected_distribution,
        device_count=device_count,
        device_name=device_name,
        compute_capability=capability,
        allowed_gpu_contracts=list(marker["allowed_gpu_contracts"]),
        observed_cuda_major=cuda_major,
        expected_cuda_major=int(marker["cuda_major"]),
    )
    native = rsc.tl.rank_genes_groups
    return (
        cp,
        rsc,
        {
            "device_count": device_count,
            "device_name": device_name,
            "compute_capability": capability,
            "validated_gpu_contract": validated,
            "rapids_singlecell": observed_version,
            "rapids_distribution": observed_distribution,
            "cuda_major": cuda_major,
            "cuda_runtime_version": runtime_version,
            "cuda_driver_version": int(cp.cuda.runtime.driverGetVersion()),
            "native_api": f"{native.__module__}.{native.__name__}",
            "native_api_signature": str(inspect.signature(native)),
            "hostname": socket.gethostname(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        },
    )
