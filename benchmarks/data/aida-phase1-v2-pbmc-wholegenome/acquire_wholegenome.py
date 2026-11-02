"""Re-read the v5 benchmark cells from Census with every gene, for protocol v6.

Purpose: give protocol v6 the same 17,394 AIDA cells, cluster labels, and embedding as
         v5, measured on the whole transcriptome instead of a 2,000-gene panel.
Run:     python benchmarks/data/aida-phase1-v2-pbmc-wholegenome/acquire_wholegenome.py OUT_H5AD
         (needs cellxgene-census and network access to the public Census bucket)
Inputs:
    ../aida-phase1-v2-pbmc-hvg2000/aida_pbmc_clustered.h5ad   pinned v5 input (sha256
        440ba41a...): obs soma_joinid, bench_cluster, annotation_scope; obsm X_pca_bench
    CELLxGENE Census 2025-11-08, Homo sapiens, the same soma_joinids
Outputs:
    OUT_H5AD   same cells in the same order; layers counts and lognorm (log1p of counts
               per 10,000 over all genes); obs and obsm copied from the pinned input
    printed    cells, genes, sha256 of OUT_H5AD
"""

import hashlib
import sys
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
PINNED = HERE.parent / "aida-phase1-v2-pbmc-hvg2000" / "aida_pbmc_clustered.h5ad"
OUT = Path(sys.argv[1])
CENSUS_VERSION = "2025-11-08"
CHUNK = 2_000  # cells per Census read; bounds peak memory

if OUT.exists():
    raise ValueError(f"refusing to overwrite {OUT}")

pinned = ad.read_h5ad(PINNED)
ids = pinned.obs["soma_joinid"].astype(np.int64).to_numpy()
if len(np.unique(ids)) != len(ids):
    raise ValueError("pinned soma_joinids are not unique")

# Read every gene for the pinned cells ---------------------------------------------------

blocks = []
with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    for start in range(0, len(ids), CHUNK):
        block = cellxgene_census.get_anndata(
            census,
            "Homo sapiens",
            obs_coords=np.sort(ids[start : start + CHUNK]),
            obs_column_names=["soma_joinid"],
            var_column_names=["feature_id", "feature_name"],
        )
        blocks.append(block)
        print(f"read {start + block.n_obs} of {len(ids)} cells")
fetched = ad.concat(blocks, merge="same")

# Align to the pinned order by soma_joinid, never by position ------------------------------
# The v10 retraction came from assigning names positionally to rows Census returned
# in joinid order. Here every row is placed by its own joinid and the result is checked.

position = {joinid: row for row, joinid in enumerate(fetched.obs["soma_joinid"].astype(np.int64))}
order = np.array([position[joinid] for joinid in ids])
fetched = fetched[order].copy()
if not np.array_equal(fetched.obs["soma_joinid"].astype(np.int64).to_numpy(), ids):
    raise ValueError("fetched cells are not in the pinned order")

counts = sp.csr_matrix(fetched.X)
totals = np.asarray(counts.sum(axis=1)).ravel()
totals[totals == 0] = 1.0
out = ad.AnnData(
    X=counts.copy(),
    obs=pinned.obs.copy(),
    var=fetched.var.set_index(fetched.var["feature_id"].astype(str)),
    obsm={key: pinned.obsm[key] for key in pinned.obsm},
    layers={
        "counts": counts,
        "lognorm": sp.csr_matrix(np.log1p(counts.multiply((1e4 / totals)[:, None]).tocsr())),
    },
)
out.var["feature_name"] = out.var["feature_name"].astype(str)
out.uns["census_version"] = CENSUS_VERSION
out.write_h5ad(OUT)

digest = hashlib.sha256(OUT.read_bytes()).hexdigest()
duplicated = int(out.var["feature_name"].duplicated().sum())
print(f"{out.n_obs} cells x {out.n_vars} genes; duplicated symbols: {duplicated}")
print(f"sha256 {digest}")
