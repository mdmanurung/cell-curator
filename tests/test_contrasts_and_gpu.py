from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from cell_curator import build_candidate_contrasts, rank_markers_gpu, validate_gpu_receipts
from cell_curator.contracts import ContractError, sha256


def _inputs():
    cells = pd.DataFrame(
        {
            "barcode": ["a", "b", "c", "d", "e", "f"],
            "scope": ["s"] * 6,
            "parent_id": ["0", "0", "0", "0", "1", "1"],
        }
    )
    registry = pd.DataFrame(
        {
            "scope": ["s"],
            "candidate_key": ["cand"],
            "parent_id": ["0"],
            "closest_sibling": ["1"],
        }
    )
    membership = pd.DataFrame(
        {
            "scope": ["s"] * 4,
            "candidate_key": ["cand"] * 4,
            "parent_id": ["0"] * 4,
            "barcode": ["a", "b", "c", "d"],
            "candidate_membership": [True, True, False, False],
        }
    )
    return cells, registry, membership


def _add_atac_plugin_lane(config: dict) -> None:
    config["input"]["assays"]["chromatin"] = {
        "modality": "atac",
        "kind": "atac",
        "feature_space": "peak",
        "capabilities": ["state"],
        "representations": {
            "tfidf": {
                "source": "X",
                "ranking_backend": "project_atac_da",
                "ontology_compatible": False,
                "detection": None,
            }
        },
    }
    config["evidence"]["lanes"].append(
        {
            "assay": "chromatin",
            "representation": "tfidf",
            "required_for_acceptance": True,
        }
    )


def test_arbitrary_lanes_expand_over_both_binary_comparators(
    tmp_path: Path, config: dict
) -> None:
    _add_atac_plugin_lane(config)
    cells, registry, membership = _inputs()
    lanes = build_candidate_contrasts.build_contrasts(
        config, cells, registry, membership, tmp_path / "contrasts"
    )
    assert len(lanes) == 4
    assert set(lanes["comparator_role"]) == {
        "parent_remainder",
        "closest_sibling",
    }
    assert set(lanes["evidence_lane_key"]) == {
        "transcriptome__logcounts",
        "chromatin__tfidf",
    }
    assert set(lanes["backend"]) == {"rapids_singlecell", "project_atac_da"}
    for path in lanes["groups_path"].unique():
        groups = pd.read_csv(path, sep="\t", dtype=str)
        assert set(groups["analysis_group"]) == {"0", "1"}


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


def _write_receipt(root: Path, lane: pd.Series) -> None:
    lane_dir = root / lane["lane_id"]
    lane_dir.mkdir(parents=True)
    wilcoxon = lane_dir / "wilcoxon.tsv.gz"
    logreg = lane_dir / "logreg.tsv.gz"
    wilcoxon.write_bytes(b"w")
    logreg.write_bytes(b"l")
    receipt = {
        "artifact_kind": "native_rapids_binary_marker_lane",
        "status": "complete",
        "lane_id": lane["lane_id"],
        "candidate_key": lane["candidate_key"],
        "comparator_role": lane["comparator_role"],
        "assay": lane["assay"],
        "assay_kind": lane["assay_kind"],
        "feature_space": lane["feature_space"],
        "representation": lane["representation"],
        "methods": ["wilcoxon", "logreg"],
        "cpu_fallback_available": False,
        "cpu_fallback_used": False,
        "runtime": {
            "device_count": 1,
            "rapids_singlecell": "9.9.9",
            "rapids_distribution": "rapids-singlecell-cu99",
            "cuda_major": 99,
            "validated_gpu_contract": {
                "name_contains": "Fixture GPU",
                "compute_capability": "9.9",
            },
        },
        "outputs": {
            "wilcoxon": {"path": str(wilcoxon), "sha256": sha256(wilcoxon)},
            "logreg": {"path": str(logreg), "sha256": sha256(logreg)},
        },
    }
    (lane_dir / "receipt.json").write_text(json.dumps(receipt))


def test_receipts_required_only_for_rapids_lanes(tmp_path: Path, config: dict) -> None:
    _add_atac_plugin_lane(config)
    cells, registry, membership = _inputs()
    lanes = build_candidate_contrasts.build_contrasts(
        config, cells, registry, membership, tmp_path / "contrasts"
    )
    root = tmp_path / "gpu"
    rapids = lanes[lanes["backend"].eq("rapids_singlecell")]
    for _, lane in rapids.iterrows():
        _write_receipt(root, lane)

    result = validate_gpu_receipts.validate_all(config, lanes, root)
    assert result["n_configured_lanes"] == 4
    assert result["n_rapids_lanes"] == 2
    assert result["rapids_lanes_per_candidate"] == 2
    assert result["cpu_fallback_used"] is False
    with pytest.raises(ContractError, match="configured binary evidence lanes"):
        validate_gpu_receipts.validate_all(config, lanes.iloc[:-1], root)
