"""The v8 comparator-to-consensus converter adds ids, flips uncertainty and keeps dropped cells."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks/scripts/v8_tools_for_consensus.py"


def test_conversion_adds_ids_confidence_and_dropped_cells(tmp_path: Path) -> None:
    query = tmp_path / "query.h5ad"
    ad.AnnData(np.zeros((4, 2)), obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"])).write_h5ad(
        query
    )
    predictions = tmp_path / "pred.tsv"
    pd.DataFrame(
        {
            "barcode": ["c1", "c2", "c3"],
            "predicted_label": ["B cell", "Unknown", "T cell"],
            "uncertainty": [0.25, 0.9, 0.0],
        }
    ).to_csv(predictions, sep="\t", index=False)
    labels = tmp_path / "label_map.tsv"
    pd.DataFrame({"label": ["B cell", "T cell"], "cl_id": ["CL:1", "CL:2"]}).to_csv(
        labels, sep="\t", index=False
    )
    out = tmp_path / "tool.tsv"
    subprocess.run(
        [sys.executable, str(SCRIPT), str(predictions), str(labels), str(query), str(out)],
        check=True,
    )
    table = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str).set_index("barcode")
    assert list(table.columns) == ["predicted_label", "predicted_cl_id", "confidence"]
    assert table.loc["c1"].tolist() == ["B cell", "CL:1", "0.75"]
    assert table.loc["c2"].tolist() == ["Unknown", "", "0.09999999999999998"]
    assert table.loc["c4"].tolist() == ["Unknown", "", "0.0"]


def test_an_unmapped_label_and_an_existing_output_both_stop_the_script(tmp_path: Path) -> None:
    query = tmp_path / "query.h5ad"
    ad.AnnData(np.zeros((1, 1)), obs=pd.DataFrame(index=["c1"])).write_h5ad(query)
    predictions = tmp_path / "pred.tsv"
    pd.DataFrame(
        {"barcode": ["c1"], "predicted_label": ["Mystery"], "uncertainty": [0.0]}
    ).to_csv(predictions, sep="\t", index=False)
    labels = tmp_path / "label_map.tsv"
    pd.DataFrame({"label": ["B cell"], "cl_id": ["CL:1"]}).to_csv(labels, sep="\t", index=False)
    out = tmp_path / "tool.tsv"

    def run() -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(predictions), str(labels), str(query), str(out)],
            capture_output=True,
            text=True,
        )

    assert "missing from the label map" in run().stderr and not out.exists()
    out.write_text("keep")
    assert "exists" in run().stderr and out.read_text() == "keep"
