"""Time the CPU marker methods on whole clusters of the AIDA whole-transcriptome slice.

Purpose: tracker row 19. How long do Wilcoxon and balanced L2 logistic regression take
on 61,497 genes, and how does it grow with the number of clusters? Nothing here changes a
call; it only measures runtime.
Run:     sbatch benchmarks/scripts/time_marker_methods.sbatch RUN_DIR
Inputs:  benchmarks/data/aida-phase1-v2-pbmc-wholegenome/aida_pbmc_wholegenome_clustered.h5ad
Outputs: RUN_DIR/timing.tsv  one row per (clusters, methods): seconds, cells, genes, fit iterations
         RUN_DIR/versions.txt
"""

import sys
import time
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from cell_curator.evidence import compute_bottom_up_markers

ROOT = Path(__file__).resolve().parents[2]
H5AD = (
    ROOT
    / "benchmarks/data/aida-phase1-v2-pbmc-wholegenome/aida_pbmc_wholegenome_clustered.h5ad"
)
LAYER = "lognorm"  # the v6 ranking layer
CLUSTER_COUNTS = [2, 4, 8]  # whole clusters, so every comparison has all of its cells
METHOD_SETS = {"wilcoxon": ("wilcoxon",), "both": ("wilcoxon", "logreg")}


def main(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=False)  # a run dir is never overwritten
    adata = ad.read_h5ad(H5AD)
    clusters = sorted(adata.obs["bench_cluster"].astype(str).unique(), key=int)
    rows = []
    for k in CLUSTER_COUNTS:
        keep = adata.obs["bench_cluster"].astype(str).isin(clusters[:k]).to_numpy()
        subset = adata[keep]
        matrix = subset.layers[LAYER]
        groups = subset.obs["bench_cluster"].astype(str)
        for name, methods in METHOD_SETS.items():
            start = time.perf_counter()
            result = compute_bottom_up_markers(
                matrix,
                groups,
                subset.var_names,
                detection=subset.layers["counts"],
                methods=methods,
            )
            seconds = time.perf_counter() - start
            fits = result.attrs["runtime"].get("logreg", {}).get("fits", {})
            rows.append(
                {
                    "clusters": k,
                    "methods": name,
                    "cells": subset.n_obs,
                    "genes": subset.n_vars,
                    "seconds": round(seconds, 1),
                    "logreg_iterations": sum(fit["iterations"] for fit in fits.values()),
                }
            )
            pd.DataFrame(rows).to_csv(
                run_dir / "timing.tsv", sep="\t", index=False
            )  # keep partial results
    (run_dir / "versions.txt").write_text(f"scanpy {sc.__version__}\n")
    sc.logging.print_versions()


if __name__ == "__main__":
    main(Path(sys.argv[1]))
