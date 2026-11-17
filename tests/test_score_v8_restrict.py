"""score_v8 --restrict drops non-listed truth cells; an unresolvable truth term is reported, not fatal."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks/scripts/score_v8.py"


def _score(tmp_path: Path, truth_ids: list[str], restrict: list[str] | None) -> dict:
    barcodes = ["c1", "c2", "c3", "c4"]
    pd.DataFrame(
        {
            "barcode": barcodes,
            "predicted_cl_id": ["CL:0000236", "CL:0000236", "CL:0000084", "CL:0000084"],
            "predicted_label": ["B cell", "B cell", "T cell", "T cell"],
            "confidence": [1.0] * 4,
        }
    ).to_csv(tmp_path / "pred.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "barcode": barcodes,
            "cell_type": ["x"] * 4,
            "cell_type_ontology_term_id": truth_ids,
        }
    ).to_csv(tmp_path / "truth.tsv", sep="\t", index=False)
    extra: list[str] = []
    if restrict is not None:
        (tmp_path / "keep.txt").write_text("\n".join(restrict) + "\n")
        extra = ["--restrict", str(tmp_path / "keep.txt")]
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(tmp_path / "pred.tsv"),
            str(tmp_path / "truth.tsv"),
            str(tmp_path / "out"),
            "--method",
            "m",
            "--dataset",
            "d",
            *extra,
        ],
        check=True,
    )
    return json.loads((tmp_path / "out/d.m.json").read_text())


def test_restrict_leaves_out_other_cells(tmp_path: Path) -> None:
    ids = ["CL:0000236", "CL:0000236", "CL:0000084", "CL:0000236"]
    score = _score(tmp_path, ids, ["c1", "c2", "c3"])
    assert score["n_truth_cells_excluded_by_restrict"] == 1
    assert score["descendant_accuracy"] == 1.0


def test_unresolvable_truth_term_is_reported(tmp_path: Path) -> None:
    ids = ["CL:0000236", "CL:0000236", "CL:0000084", "CL:9999999"]
    score = _score(tmp_path, ids, None)
    assert score["n_truth_cells_unmapped"] == 1
    assert score["unmapped_truth_cl_ids"] == ["CL:9999999"]
