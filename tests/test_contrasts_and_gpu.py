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
            "ranking_options": [
                "use_continuity=False applied natively",
                "multi_gpu=False applied natively",
            ],
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


class _FakeTools:
    """Stand-in for rapids_singlecell.tl with a controllable signature."""

    def __init__(self, *, with_wilcoxon: bool, extra: tuple[str, ...] = ()) -> None:
        if not with_wilcoxon:
            return

        if extra == ("use_continuity", "multi_gpu"):

            def rank_genes_groups(  # type: ignore[no-redef]
                adata, groupby, *, tie_correct=False, use_continuity=False,
                multi_gpu=None, **kwds,
            ): ...
        else:

            def rank_genes_groups(adata, groupby, *, tie_correct=False, **kwds): ...

        self.rank_genes_groups = rank_genes_groups


class _FakeRsc:
    def __init__(self, tools: _FakeTools) -> None:
        self.tl = tools


def test_a_release_without_the_ranking_entry_point_names_itself() -> None:
    """Better than a bare AttributeError from deep inside the lane."""

    rsc = _FakeRsc(_FakeTools(with_wilcoxon=False))
    with pytest.raises(
        rank_markers_gpu.GpuContractError, match="rank_genes_groups is required"
    ):
        rank_markers_gpu.native_ranking_callable(rsc)


def test_only_declared_parameters_count_as_supported() -> None:
    """**kwds would swallow an unknown option without applying it."""

    older = _FakeRsc(_FakeTools(with_wilcoxon=True))
    supported = rank_markers_gpu.supported_ranking_options(older.tl.rank_genes_groups)
    assert "tie_correct" in supported
    assert "use_continuity" not in supported

    newer = _FakeRsc(_FakeTools(with_wilcoxon=True, extra=("use_continuity", "multi_gpu")))
    assert {"use_continuity", "multi_gpu"}.issubset(
        rank_markers_gpu.supported_ranking_options(newer.tl.rank_genes_groups)
    )


def test_an_option_requested_off_is_satisfied_by_its_absence() -> None:
    include, note = rank_markers_gpu.resolve_optional_ranking_option(
        name="use_continuity",
        requested=False,
        supported=frozenset({"tie_correct"}),
        installed_version="0.14.1",
    )
    assert include is False
    assert "requested off" in note and "0.14.1" in note


def test_an_option_requested_on_fails_closed_when_it_cannot_be_applied() -> None:
    """Silently dropping it would produce statistics nobody asked for."""

    with pytest.raises(rank_markers_gpu.GpuContractError, match="cannot be applied"):
        rank_markers_gpu.resolve_optional_ranking_option(
            name="use_continuity",
            requested=True,
            supported=frozenset({"tie_correct"}),
            installed_version="0.14.1",
        )


def test_a_supported_option_is_passed_through(config: dict) -> None:
    include, note = rank_markers_gpu.resolve_optional_ranking_option(
        name="use_continuity",
        requested=True,
        supported=frozenset({"use_continuity"}),
        installed_version="0.16.1",
    )
    assert include is True
    assert "applied natively" in note
    # The configuration's own request feeds the same resolution path.
    assert dict(rank_markers_gpu.optional_ranking_requests(config["markers"])) == {
        "use_continuity": False,
        "multi_gpu": False,
    }
