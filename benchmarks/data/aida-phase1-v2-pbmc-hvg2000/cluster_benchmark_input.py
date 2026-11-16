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

import math
import sys
from pathlib import Path

import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.special import comb


# Copied from cell_curator.refinement, removed with the package's own clustering.
def _knn_graph(embedding: np.ndarray, *, k: int) -> nx.Graph:
    if embedding.ndim != 2 or len(embedding) < 2:
        raise ValueError("the embedding needs at least two cells")
    effective_k = min(k + 1, len(embedding))
    distances, indices = cKDTree(embedding).query(embedding, k=effective_k)
    if indices.ndim == 1:
        indices = indices[:, None]
        distances = distances[:, None]
    graph = nx.Graph()
    graph.add_nodes_from(range(len(embedding)))
    positive = distances[:, 1:][distances[:, 1:] > 0]
    scale = float(np.median(positive)) if positive.size else 1.0
    for source in range(len(embedding)):
        for distance, target in zip(distances[source, 1:], indices[source, 1:], strict=True):
            if source == int(target):
                continue
            weight = math.exp(-float(distance) / max(scale, 1e-12))
            if graph.has_edge(source, int(target)):
                graph[source][int(target)]["weight"] = max(
                    graph[source][int(target)]["weight"], weight
                )
            else:
                graph.add_edge(source, int(target), weight=weight)
    return graph


def graph_partition(
    embedding: np.ndarray,
    *,
    resolution: float,
    seed: int,
    k: int,
) -> np.ndarray:
    graph = _knn_graph(embedding, k=k)
    communities = nx.community.louvain_communities(
        graph,
        weight="weight",
        resolution=float(resolution),
        seed=int(seed),
    )
    ordered = sorted(
        (sorted(group) for group in communities), key=lambda group: (group[0], len(group))
    )
    labels = np.empty(len(embedding), dtype=int)
    for label, members in enumerate(ordered):
        labels[members] = label
    return labels


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
