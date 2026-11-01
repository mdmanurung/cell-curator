"""Re-derive the benchmark clustering from the recipe in benchmarks/README.md.

Purpose: turn the prose recipe for `bench_cluster` into code and measure how
         closely it reproduces the pinned clustered input, which no surviving
         script produced.
Run:     uv run python benchmarks/data/aida-phase1-v2-pbmc-hvg2000/cluster_benchmark_input.py OUT_TSV
Inputs:
    aida_pbmc_clustered.h5ad   pinned v3-v5 input (sha256 440ba41a...); layers lognorm,
                               obs bench_cluster [benchmarks/README.md, recipe only]
Outputs:
    OUT_TSV                    barcode, recomputed cluster, pinned cluster
    printed                    cluster counts and adjusted Rand index against the pinned
                               clusters
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.special import comb

from cell_curator.refinement import graph_partition

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1])
# The recipe in benchmarks/README.md, "The input caps this evaluation".
N_PCS = 50
N_PCS_GRAPH = 30
K_NEIGHBOURS = 15
RESOLUTION = 1.0
SEED = 0

if OUT.exists():
    raise ValueError(f"refusing to overwrite {OUT}")

# Scale and project ---------------------------------------------------------------
# Unit-variance scaling of the log-normalized layer, then PCA by SVD. The README does
# not name the layer or the PCA implementation; these are the choices made here.

adata = ad.read_h5ad(HERE / "aida_pbmc_clustered.h5ad")
values = adata.layers["lognorm"]
values = values.toarray() if hasattr(values, "toarray") else np.asarray(values)
values = values.astype(np.float64)
spread = values.std(axis=0)
scaled = (values - values.mean(axis=0)) / np.where(spread > 0, spread, 1.0)
left, singular, _ = np.linalg.svd(scaled, full_matrices=False)
pcs = left[:, :N_PCS] * singular[:N_PCS]

# Cluster -------------------------------------------------------------------------

recomputed = graph_partition(
    pcs[:, :N_PCS_GRAPH], resolution=RESOLUTION, seed=SEED, k=K_NEIGHBOURS
)
table = pd.DataFrame(
    {
        "barcode": adata.obs_names,
        "recomputed": recomputed.astype(str),
        "pinned": adata.obs["bench_cluster"].astype(str).to_numpy(),
    }
)
table.to_csv(OUT, sep="\t", index=False)

# Compare with the pinned clusters -------------------------------------------------

contingency = pd.crosstab(table["recomputed"], table["pinned"]).to_numpy()
pairs = comb(contingency, 2).sum()
rows = comb(contingency.sum(axis=1), 2).sum()
columns = comb(contingency.sum(axis=0), 2).sum()
expected = rows * columns / comb(len(table), 2)
ari = (pairs - expected) / ((rows + columns) / 2 - expected)
print(f"recomputed clusters: {table['recomputed'].nunique()}")
print(f"pinned clusters:     {table['pinned'].nunique()}")
print(f"adjusted Rand index: {ari:.4f}")
