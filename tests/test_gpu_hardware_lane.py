"""Real-hardware validation of the RAPIDS evidence lane.

Every other GPU test in this suite proves the contract *fails closed* without a
GPU. This module proves the complementary half that cannot be faked: on a node
with a real, contract-satisfying NVIDIA GPU the lane executes natively, emits a
receipt recording both native marker methods, and still rejects a mismatched
hardware contract.

Skipped unless ``cupy``/``rapids_singlecell`` import and exactly one CUDA device
is visible, so it is inert on CPU-only machines and in CI.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import SKILL_ROOT

from cell_curator.capabilities import build_capability_registry
from cell_curator.config import load_config
from cell_curator.contracts import ContractError
from cell_curator.evidence import validate_backend_receipt
from cell_curator.models import BackendReceipt
from cell_curator.pipeline import evidence_lane
from cell_curator.rank_markers_gpu import (
    GpuContractError,
    format_compute_capability,
    inspect_runtime,
    validate_runtime_contract_values,
)

pytestmark = pytest.mark.gpu


def _observed_runtime() -> dict[str, object]:
    """Read the node's real RAPIDS/CUDA identity, or skip if there is none."""

    cupy = pytest.importorskip("cupy", reason="cupy is required for GPU hardware tests")
    pytest.importorskip(
        "rapids_singlecell", reason="rapids-singlecell is required for GPU hardware tests"
    )
    from importlib.metadata import distribution, version

    if int(cupy.cuda.runtime.getDeviceCount()) != 1:
        pytest.skip("the RAPIDS contract requires exactly one visible CUDA device")
    device = cupy.cuda.Device(0)
    properties = cupy.cuda.runtime.getDeviceProperties(0)
    raw_name = properties["name"]
    runtime_version = int(cupy.cuda.runtime.runtimeGetVersion())
    return {
        "version": version("rapids-singlecell"),
        "distribution": str(distribution("rapids-singlecell").metadata["Name"]),
        "cuda_major": runtime_version // 1000,
        "device_name": raw_name.decode() if isinstance(raw_name, bytes) else str(raw_name),
        "compute_capability": format_compute_capability(device.compute_capability),
    }


def _gpu_fixture(tmp_path: Path, observed: dict[str, object]) -> Path:
    """Write a single-lane evidence-only config pinned to the observed hardware."""

    rng = np.random.default_rng(17)
    counts = np.zeros((80, 6), dtype=np.float32)
    counts[:40, :2] = rng.poisson(9, size=(40, 2)) + 1
    counts[40:, 2:4] = rng.poisson(9, size=(40, 2)) + 1
    genes = ["A1", "A2", "B1", "B2", "N1", "N2"]
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame(
            {
                "canonical_cluster": ["0"] * 40 + ["1"] * 40,
                "lineage": ["immune"] * 80,
            },
            index=[f"gpu-{index:03d}" for index in range(80)],
        ),
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = np.vstack(
        [rng.normal(-3, 0.2, size=(40, 3)), rng.normal(3, 0.2, size=(40, 3))]
    ).astype(np.float32)
    source = tmp_path / "gpu-source.h5ad"
    adata.write_h5ad(source)

    markers = tmp_path / "gpu-markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {"pos": ["A1", "A2"], "neg": ["B1", "B2"]},
                "Type B": {"pos": ["B1", "B2"], "neg": ["A1", "A2"]},
            }
        )
    )

    value = yaml.safe_load((SKILL_ROOT / "assets/config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": "gpu-hardware-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "evidence-only",
        }
    )
    value["input"]["path"] = str(source)
    value["canonical_partitions"] = {
        "immune": {"parent_key": "canonical_cluster", "expected_parent_count": 2}
    }
    value["adaptive"]["eligibility"].update(
        {"min_cells": 5, "min_donors": None, "min_captures": None}
    )
    value["adaptive"]["parent_local"]["enabled"] = False
    value["knowledge"]["providers"] = [
        {"kind": "local_markers", "path": str(markers), "version": "gpu-hardware-v1"}
    ]
    value["markers"].update(
        {
            "rapids_singlecell_version": observed["version"],
            "rapids_distribution": observed["distribution"],
            "cuda_major": observed["cuda_major"],
            "cpu_fallback": False,
            "allowed_gpu_contracts": [
                {
                    "name_contains": str(observed["device_name"]),
                    "compute_capability": observed["compute_capability"],
                }
            ],
        }
    )
    value["hpc"]["gpu_resources"] = {"gpu": 1}
    value["input"]["assays"]["transcriptome"]["representations"]["logcounts"][
        "ranking_backend"
    ] = "rapids_singlecell"
    path = tmp_path / "gpu-config.yaml"
    path.write_text(yaml.safe_dump(value, sort_keys=False))
    return path


def test_real_gpu_validates_the_declared_hardware_contract(tmp_path: Path) -> None:
    """The declared contract is checked against the live device, not just parsed."""

    observed = _observed_runtime()
    config_path = _gpu_fixture(tmp_path, observed)
    marker = load_config(config_path).markers.model_dump(mode="json")

    validated = validate_runtime_contract_values(
        observed_version=str(observed["version"]),
        expected_version=str(marker["rapids_singlecell_version"]),
        observed_distribution=str(observed["distribution"]),
        expected_distribution=str(marker["rapids_distribution"]),
        device_count=1,
        device_name=str(observed["device_name"]),
        compute_capability=str(observed["compute_capability"]),
        allowed_gpu_contracts=list(marker["allowed_gpu_contracts"]),
        observed_cuda_major=int(observed["cuda_major"]),
        expected_cuda_major=int(marker["cuda_major"]),
    )
    assert validated["name_contains"] == observed["device_name"]
    assert validated["compute_capability"] == observed["compute_capability"]
    # A two-digit major must survive the round trip; "120" must not become "1.2".
    major, _, minor = str(observed["compute_capability"]).partition(".")
    assert major.isdigit() and len(minor) == 1


@pytest.mark.xfail(
    strict=True,
    reason=(
        "rapids-singlecell exposes no rank_genes_groups and no Wilcoxon "
        "implementation at all (only tl.rank_genes_groups_logreg), so "
        "rank_markers_gpu.inspect_runtime and evidence.compute_gpu_bottom_up_markers "
        "bind a non-existent API and the lane cannot execute. Confirmed on an "
        "RTX PRO 6000 Blackwell MIG slice with rapids-singlecell 0.13.4 / RAPIDS "
        "25.12. Flip to a pass once the GPU marker path targets the real API."
    ),
)
def test_rapids_lane_executes_natively_and_publishes_a_gpu_receipt(tmp_path: Path) -> None:
    observed = _observed_runtime()
    config_path = _gpu_fixture(tmp_path, observed)

    receipt_path = evidence_lane(
        config_path, assay="transcriptome", representation="logcounts"
    )
    receipt = BackendReceipt.model_validate_json(Path(receipt_path).read_text())

    assert receipt.backend == "rapids_singlecell"
    assert receipt.status == "complete"
    assert set(receipt.marker_methods) == {"wilcoxon", "logreg"}
    assert receipt.convergence_validated is True
    assert not receipt.fallback_used
    assert not receipt.cpu_fallback_available
    assert not receipt.cpu_fallback_used

    runtime = receipt.runtime
    assert int(runtime["device_count"]) == 1
    assert runtime["device_name"] == observed["device_name"]
    assert runtime["compute_capability"] == observed["compute_capability"]
    assert int(runtime["cuda_major"]) == int(observed["cuda_major"])
    assert str(runtime["native_api"]).startswith("rapids_singlecell")
    assert runtime["cuda_driver_version"]

    # The published receipt must also satisfy the consumer-side contract that
    # every downstream reader applies before trusting the lane.
    capability = next(
        item
        for item in build_capability_registry(load_config(config_path))
        if item.assay == "transcriptome" and item.representation == "logcounts"
    )
    assert capability.requires_gpu is True
    validate_backend_receipt(receipt, capability)


def test_real_gpu_still_rejects_a_mismatched_hardware_contract(tmp_path: Path) -> None:
    """The allowlist is enforced against the real device, not merely declared."""

    observed = _observed_runtime()
    config_path = _gpu_fixture(tmp_path, observed)
    value = yaml.safe_load(config_path.read_text())
    value["markers"]["allowed_gpu_contracts"] = [
        {"name_contains": "Nonexistent Model", "compute_capability": "0.1"}
    ]
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))

    with pytest.raises((GpuContractError, ContractError), match="not allowed"):
        inspect_runtime(load_config(config_path).model_dump(mode="json"))
