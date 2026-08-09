from __future__ import annotations

import pytest

from cell_curator import rank_markers_gpu
from cell_curator.contracts import ContractError
from cell_curator.evidence import validate_backend_receipt
from cell_curator.models import (
    AssayKind,
    BackendReceipt,
    Capability,
    CapabilityLevel,
)


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"device_count": 0}, "exactly one"),
        ({"device_name": "CPU"}, "not allowed"),
        ({"observed_version": "0.0.0"}, "version mismatch"),
    ],
)
def test_gpu_contract_fails_closed(config: dict, updates: dict, message: str) -> None:
    values = {
        "observed_version": "9.9.9",
        "expected_version": "9.9.9",
        "observed_distribution": "rapids-singlecell-cu99",
        "expected_distribution": "rapids-singlecell-cu99",
        "device_count": 1,
        "device_name": "NVIDIA Fixture GPU",
        "compute_capability": "9.9",
        "allowed_gpu_contracts": config["markers"]["allowed_gpu_contracts"],
        "observed_cuda_major": 99,
        "expected_cuda_major": 99,
    }
    values.update(updates)
    with pytest.raises(ContractError, match=message):
        rank_markers_gpu.validate_runtime_contract_values(**values)


def _receipt(**updates) -> BackendReceipt:
    value = {
        "run_id": "gpu-fixture",
        "lane_key": "rna__logcounts",
        "backend": "rapids_singlecell",
        "capability_level": "annotation-capable",
        "status": "complete",
        "input_sha256": "a" * 64,
        "config_sha256": "b" * 64,
        "marker_methods": ["wilcoxon", "logreg"],
        "cpu_fallback_available": False,
        "cpu_fallback_used": False,
        "fallback_used": False,
        "convergence_validated": True,
        "runtime": {
            "device_count": 1,
            "device_name": "NVIDIA Fixture GPU",
            "compute_capability": "9.9",
            "validated_gpu_contract": {
                "name_contains": "Fixture GPU",
                "compute_capability": "9.9",
            },
            "rapids_singlecell": "9.9.9",
            "rapids_distribution": "rapids-singlecell-cu99",
            "cuda_major": 99,
            "cuda_runtime_version": 99000,
            "cuda_driver_version": 99000,
            "native_api": "rapids_singlecell.tl.rank_genes_groups",
        },
    }
    value.update(updates)
    return BackendReceipt.model_validate(value)


def test_integrated_gpu_receipt_requires_both_native_methods_and_no_fallback() -> None:
    capability = Capability(
        assay="rna",
        representation="logcounts",
        assay_kind=AssayKind.RNA,
        feature_space="gene",
        level=CapabilityLevel.ANNOTATION,
        claims=("cell_identity",),
        backend="rapids_singlecell",
        requires_gpu=True,
    )
    validate_backend_receipt(_receipt(), capability)
    with pytest.raises(ContractError, match="Wilcoxon and balanced logreg"):
        validate_backend_receipt(_receipt(marker_methods=["wilcoxon"]), capability)
    with pytest.raises(ContractError, match="prohibited CPU fallback"):
        validate_backend_receipt(_receipt(cpu_fallback_used=True), capability)


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("86", "8.6"),
        ("75", "7.5"),
        ("90", "9.0"),
        # Blackwell sm_120: positional splitting would yield "1.2" and silently
        # compare against allowed_gpu_contracts as a different device.
        ("120", "12.0"),
        ("100", "10.0"),
        ((12, 0), "12.0"),
        ("12.0", "12.0"),
    ],
)
def test_compute_capability_formats_two_digit_majors(raw: object, expected: str) -> None:
    assert rank_markers_gpu.format_compute_capability(raw) == expected


@pytest.mark.parametrize("raw", ["", "x", "9", (1, 2, 3)])
def test_compute_capability_rejects_unreadable_values(raw: object) -> None:
    with pytest.raises(rank_markers_gpu.GpuContractError, match="unreadable compute capability"):
        rank_markers_gpu.format_compute_capability(raw)
