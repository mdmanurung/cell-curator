# Purpose: per-cell CellTypist labels (no majority voting) for one public model on the PBMC3K
#          input, mapped to Cell Ontology ids with the shared mapping table.
# Run:     /exports/para-lipg-hpc/mdmanurung/cc-v8-work/envs/public-models/venv/bin/python \
#              02a_run-celltypist.py MODEL RUN_DIR
# Inputs:  RUN_DIR/pbmc3k.h5ad                      layers logcounts (log1p CP10k) [01_prep-pbmc3k.py]
#          ~/.celltypist/data/models/MODEL.pkl      public model
#          public_models/celltypist/mapping.tsv     model, celltypist_label, cl_id, cl_name
# Outputs: RUN_DIR/tools/celltypist.tsv   barcode, predicted_label, predicted_cl_id, confidence
#          RUN_DIR/tools/celltypist.raw.tsv  barcode, celltypist_label, confidence

import sys
from pathlib import Path

import anndata as ad
import celltypist
import pandas as pd

MAPPING = Path(
    "/exports/para-lipg-hpc/mdmanurung/cc-v8-work/public_models/celltypist/mapping.tsv"
)

model, run_dir = sys.argv[1], Path(sys.argv[2])
out = run_dir / "tools"
out.mkdir(exist_ok=True)

# CellTypist, per cell -------------------------------------------------------------
adata = ad.read_h5ad(run_dir / "pbmc3k.h5ad")
query = ad.AnnData(
    X=adata.layers["logcounts"], obs=pd.DataFrame(index=adata.obs_names), var=adata.var[[]]
)
pred = celltypist.annotate(
    query, model=model + ".pkl", majority_voting=False, mode="best match"
)
raw = pd.DataFrame(
    {
        "barcode": pred.predicted_labels.index,
        "celltypist_label": pred.predicted_labels["predicted_labels"].to_numpy(),
        "confidence": pred.probability_matrix.max(axis=1).to_numpy(),
    }
)
raw.to_csv(out / "celltypist.raw.tsv", sep="\t", index=False)

# Map labels to CL ids; a label with no mapping row must fail -----------------------
mapping = pd.read_csv(MAPPING, sep="\t", keep_default_na=False)
mapping = mapping[mapping["model"] == model].drop(columns="model")
table = raw.merge(mapping, on="celltypist_label", how="left", validate="m:1")
if table["source"].isna().any():
    raise ValueError(
        f"labels with no mapping row: {set(table.loc[table.source.isna(), 'celltypist_label'])}"
    )
table = table.rename(columns={"cl_id": "predicted_cl_id", "cl_name": "predicted_label"})
table[["barcode", "predicted_label", "predicted_cl_id", "confidence"]].to_csv(
    out / "celltypist.tsv", sep="\t", index=False
)
