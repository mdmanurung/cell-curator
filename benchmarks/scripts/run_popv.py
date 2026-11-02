"""Run popV trained on v5 donor split A and annotate split B.

Purpose: the popV comparator for protocol v5. OnClass, CellTypist, and scANVI are left
         out of popV's ensemble, so popV is not an average of comparators scored on
         their own; the seven remaining methods run.
Run:     python benchmarks/scripts/run_popv.py SPLITS_DIR OUT_DIR   (on a GPU node)
Inputs:
    SPLITS_DIR/A.h5ad, SPLITS_DIR/B.h5ad   [export_splits.py]; layers counts, obs cell_type
Outputs (OUT_DIR, must not exist):
    cell_predictions.tsv   barcode, predicted_label, uncertainty
                           uncertainty = 1 - agreeing methods / methods run
    models/                popV's trained models
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import popv

SPLITS = Path(sys.argv[1])
OUT = Path(sys.argv[2])
METHODS = [
    "KNN_BBKNN",
    "KNN_HARMONY",
    "KNN_SCANORAMA",
    "KNN_SCVI",
    "Random_Forest",
    "Support_Vector",
    "XGboost",
]
BATCH_KEY = "donor_id"  # popV corrects batch by this key in both splits

# UMAP embeddings are for popV's plots, not its predictions. Off, because scanpy 1.11
# cannot read the sparse graph bbknn 1.6 builds for that UMAP.
popv.settings.compute_umap_embedding = False

OUT.mkdir(parents=True)
reference = ad.read_h5ad(SPLITS / "A.h5ad")
query = ad.read_h5ad(SPLITS / "B.h5ad")
for part in (reference, query):
    part.X = part.layers["counts"].copy()

# Train and annotate ------------------------------------------------------------------
# hvg=None keeps the panel's 2,000 genes; popV's default of 4,000 exceeds them.

processed = popv.preprocessing.Process_Query(
    query,
    reference,
    ref_labels_key="cell_type",
    ref_batch_key=BATCH_KEY,
    query_batch_key=BATCH_KEY,
    cl_obo_folder=False,
    prediction_mode="retrain",
    save_path_trained_models=str(OUT / "models"),
    hvg=None,
).adata
popv.annotation.annotate_data(processed, methods=METHODS)

# Write ------------------------------------------------------------------------------

called = processed[processed.obs["_dataset"].eq("query")].obs
score = called["popv_prediction_score"].astype(float)
table = pd.DataFrame(
    {
        "barcode": called.index.str.removesuffix("-query"),
        "predicted_label": called["popv_prediction"].astype(str).to_numpy(),
        "uncertainty": (1 - score / len(METHODS)).to_numpy(),
    }
)
if not set(table["barcode"]).issubset(query.obs_names):
    raise ValueError("popV changed query barcodes; fix the suffix rule before scoring")
table.to_csv(OUT / "cell_predictions.tsv", sep="\t", index=False)
print(f"{len(table)} of {query.n_obs} query cells annotated")
