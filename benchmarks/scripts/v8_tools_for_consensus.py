"""Turn a comparator's predictions into a tool table for the annotation-consensus skill.

Purpose: the consensus skill reads one TSV per tool (barcode, predicted_label, predicted_cl_id,
         confidence). The v8 comparators write barcode, predicted_label, uncertainty. This adds
         the Cell Ontology id from the training label map and sets confidence = 1 - uncertainty,
         the rule that produced the dev-set tool files (checked in test_v8_tools_for_consensus.py).
Run:     .venv/bin/python benchmarks/scripts/v8_tools_for_consensus.py PREDICTIONS LABEL_MAP QUERY_H5AD OUT_TSV
Inputs:  PREDICTIONS  barcode, predicted_label, uncertainty   [v8_celltypist.py, run_popv.py, v8_singler.R]
                      CellTypist: use majority/cell_predictions.tsv (the dev-set tool file was the
                      cluster-majority lane), not per_cell
         LABEL_MAP    TRAIN_DIR/label_map.tsv: label, cl_id   [v8_prepare_training.py]
         QUERY_H5AD   input.h5ad; a cell the tool dropped (popV drops some in preprocessing) is added
                      as Unknown with a blank id and confidence 0, as in the dev-set files
Outputs: OUT_TSV      barcode, predicted_label, predicted_cl_id, confidence (must not exist)
The label Unknown (SingleR's own abstention) keeps a blank id. Any other label missing from the
label map stops the script.
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd

predictions, label_map, query, out = (Path(p) for p in sys.argv[1:5])
if out.exists():
    raise SystemExit(f"{out} exists; this script never overwrites")
table = pd.read_csv(predictions, sep="\t", keep_default_na=False, dtype={"barcode": str})
ids = pd.read_csv(label_map, sep="\t", keep_default_na=False, dtype=str).set_index("label")[
    "cl_id"
]
unknown = sorted(set(table["predicted_label"]) - set(ids.index) - {"Unknown"})
if unknown:
    raise SystemExit(f"labels missing from the label map: {unknown[:5]}")
table["predicted_cl_id"] = table["predicted_label"].map(ids).fillna("")
table["confidence"] = 1 - table["uncertainty"].astype(float)
table = table[["barcode", "predicted_label", "predicted_cl_id", "confidence"]]
missing = ad.read_h5ad(query, backed="r").obs_names.difference(table["barcode"])
dropped = pd.DataFrame(
    {"barcode": missing, "predicted_label": "Unknown", "predicted_cl_id": "", "confidence": 0.0}
)
pd.concat([table, dropped], ignore_index=True).to_csv(out, sep="\t", index=False)
