"""The real-data tutorial must run to labels, because it is what a new user tries first."""

import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("scanpy")

SCRIPT = Path(__file__).resolve().parents[1] / "examples" / "pbmc3k_tutorial.py"


def test_the_pbmc3k_tutorial_runs_to_a_label_table(tmp_path: Path) -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=1800,
    )
    assert completed.returncode == 0, completed.stderr[-4000:]
    assert "scanpy_name" in completed.stdout
    labels = tmp_path / "results" / "pbmc3k-tutorial-v1" / "labels" / "L1.tsv"
    assert len(labels.read_text().splitlines()) == 9  # header plus 8 clusters
