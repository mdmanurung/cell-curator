from __future__ import annotations

from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPO_ROOT / "skills" / "cell-curator"


@pytest.fixture
def config(tmp_path: Path) -> dict:
    value = yaml.safe_load((SKILL_ROOT / "assets/config.template.yaml").read_text())
    value["run"].update(
        {
            "run_id": "fixture-v1",
            "output_root": str(tmp_path / "results"),
            "mode": "evidence-only",
        }
    )
    value["adaptive"]["eligibility"].update(
        {"min_cells": 2, "min_donors": 2, "min_captures": 2}
    )
    value["input"]["keys"].update(
        {
            "donor": "donor_id",
            "capture": "capture_id",
            "doublet_class": "doublet_class",
        }
    )
    value["adaptive"]["audit"].update(
        {
            "donor_max_fraction": 0.6,
            "capture_max_fraction": 0.6,
            "doublet_fraction_max": 0.15,
        }
    )
    value["adaptive"]["stored_candidates"].update(
        {
            "min_parent_purity": 0.9,
            "min_parent_fraction": 0.2,
            "max_parent_fraction": 0.95,
            "adjacent_jaccard_min": 0.8,
            "min_adjacent_matches": 1,
        }
    )
    value["markers"].update(
        {
            "rapids_singlecell_version": "9.9.9",
            "rapids_distribution": "rapids-singlecell-cu99",
            "cuda_major": 99,
            "allowed_gpu_contracts": [
                {"name_contains": "Fixture GPU", "compute_capability": "9.9"}
            ],
        }
    )
    # Exercise the fail-closed GPU lane independently of the CPU-safe template.
    value["input"]["assays"]["transcriptome"]["representations"]["logcounts"][
        "ranking_backend"
    ] = "rapids_singlecell"
    return value


@pytest.fixture
def config_path(tmp_path: Path, config: dict) -> Path:
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path
