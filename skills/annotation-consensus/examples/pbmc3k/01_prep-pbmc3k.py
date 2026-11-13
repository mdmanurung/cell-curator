# Purpose: build the PBMC3K input of the annotation-consensus showcase: counts, log1p CP10k,
#          PCA and the scanpy tutorial's louvain clusters as opaque cluster codes.
# Run:     .venv/bin/python 01_prep-pbmc3k.py RUN_DIR      (RUN_DIR must not exist)
# Inputs:  sc.datasets.pbmc3k()            raw counts, 2,700 cells (downloaded)
#          sc.datasets.pbmc3k_processed()  scanpy tutorial clusters, 2,638 cells (downloaded)
# Outputs: RUN_DIR/pbmc3k.h5ad   layers counts + logcounts, obsm X_pca, obs cluster (code),
#                                donor_id and scanpy_name (kept for the comparison only, never voted on)
#          RUN_DIR/prep.log      cell and gene counts

import sys
from pathlib import Path

import scanpy as sc

# Parameters -----------------------------------------------------------------------
N_PCS = 50  # the scanpy tutorial's PCA size
MIN_CELLS = 3  # the scanpy tutorial's gene filter
SEED = 0

run_dir = Path(sys.argv[1])
run_dir.mkdir(parents=True)
sc.settings.datasetdir = run_dir / "data"

# Load counts and tutorial clusters ------------------------------------------------
# Counts come from the raw dataset; the tutorial filtered cells, so keep only those.

adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
names = sc.datasets.pbmc3k_processed().obs["louvain"]
n_raw = adata.n_obs
adata = adata[adata.obs_names.isin(names.index)].copy()
names = names.loc[adata.obs_names]
adata.obs["cluster"] = names.cat.codes.astype(str).to_numpy()  # codes: no names reach the run
adata.obs["scanpy_name"] = names.to_numpy()
adata.obs["donor_id"] = "pbmc3k"  # one donor; popV needs a batch key

# Normalize and embed --------------------------------------------------------------

n_genes = adata.n_vars
sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
sc.pp.pca(adata, n_comps=N_PCS, random_state=SEED)

adata.write_h5ad(run_dir / "pbmc3k.h5ad")
(run_dir / "prep.log").write_text(
    f"cells {adata.n_obs} of {n_raw}\ngenes {adata.n_vars} of {n_genes}\n"
    f"clusters {adata.obs['cluster'].nunique()}\n"
)
sc.logging.print_versions()
