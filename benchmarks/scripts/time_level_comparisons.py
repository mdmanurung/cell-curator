"""Time one whole annotation level on the AIDA whole-transcriptome slice.

Purpose: protocol v7 cap (Q7). The corpus workflow runs one parent-remainder comparison over
all clusters, then one pair comparison per cluster against its nearest cluster, each with
Wilcoxon and logistic regression. This replays that structure with the same function so the
wall time of a real level is measured, not extrapolated.
Run:     sbatch benchmarks/scripts/time_level_comparisons.sbatch RUN_DIR
Inputs:  benchmarks/data/aida-phase1-v2-pbmc-wholegenome/aida_pbmc_wholegenome_clustered.h5ad
Outputs: RUN_DIR/steps.tsv  one row per comparison: step, cluster, cells, seconds (written as it goes)
"""

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from cell_curator.evidence import compute_bottom_up_markers

ROOT = Path(__file__).resolve().parents[2]
H5AD = (
    ROOT
    / "benchmarks/data/aida-phase1-v2-pbmc-wholegenome/aida_pbmc_wholegenome_clustered.h5ad"
)
LAYER = "lognorm"  # the v6 ranking layer
EMBEDDING = "X_pca_bench"  # the v6 embedding, used for the nearest cluster
METHODS = ("wilcoxon", "logreg")


def main(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=False)  # a run dir is never overwritten
    adata = ad.read_h5ad(H5AD)
    groups = adata.obs["bench_cluster"].astype(str)
    clusters = sorted(groups.unique(), key=int)
    embedding = np.asarray(adata.obsm[EMBEDDING])
    centroids = {c: embedding[(groups == c).to_numpy()].mean(axis=0) for c in clusters}
    steps = []

    def timed(step: str, cluster: str, keep: np.ndarray) -> None:
        subset = adata[keep]
        start = time.perf_counter()
        compute_bottom_up_markers(
            subset.layers[LAYER],
            subset.obs["bench_cluster"].astype(str),
            subset.var_names,
            detection=subset.layers["counts"],
            top_n=None,
            methods=METHODS,
        )
        steps.append(
            {
                "step": step,
                "cluster": cluster,
                "cells": subset.n_obs,
                "seconds": round(time.perf_counter() - start, 1),
            }
        )
        pd.DataFrame(steps).to_csv(
            run_dir / "steps.tsv", sep="\t", index=False
        )  # keep partial results

    timed("parent_remainder", "all", np.ones(adata.n_obs, dtype=bool))
    for cluster in clusters:
        nearest = min(
            (c for c in clusters if c != cluster),
            key=lambda c: (float(np.linalg.norm(centroids[cluster] - centroids[c])), c),
        )
        timed(
            "closest_sibling",
            f"{cluster}-{nearest}",
            groups.isin([cluster, nearest]).to_numpy(),
        )


if __name__ == "__main__":
    main(Path(sys.argv[1]))
