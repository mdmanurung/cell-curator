"""Build the annotation-consensus demo fixture from AIDA split A, deterministically.

Purpose: draw about 600 cells (seeded, stratified so every bench_cluster is present) and write
         the subset of every input the pipeline reads. Truth goes to a separate file that only
         the score step reads.
Run:     .venv/bin/python skills/annotation-consensus/examples/make_demo.py
Inputs:
  prep-aida_split_a-2026-11-09/input.h5ad        obs bench_cluster, obsm X_pca [prep]
  c004_celllevel/cell_votes.tsv                  per-cell marker votes [c004 01_cell_votes.py]
  mvp-aida-2026-11-10/predictions_{celltypist,popv}.tsv   tool predictions [mvp-aida]
  c012_harmonizer/preds/aida_predictions_singler.tsv      SingleR converted [c012]
  runs/v8-dev-aida-split-a-003/plan/vocabulary.tsv        corpus vocabulary [run 003]
  runs/v8-dev-aida-split-a-003/annotation/levels/L1/markers.tsv   marker ranking [run 003]
  prep-aida_split_a-2026-11-09/scoring/truth.tsv  truth, written to truth.tsv for scoring only
Outputs (examples/demo/):
  cells.h5ad, cell_votes.tsv, predictions_{celltypist,popv,singler}.tsv, vocabulary.tsv,
  markers.tsv, truth.tsv
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Parameters -----------------------------------------------------------------------
SEED = 20261112  # the draw is reproducible; change it and the demo changes
TARGET_CELLS = 600  # about this many cells in total
MIN_PER_CLUSTER = 10  # every cluster is present, with at least this many cells (or all of it)
TOP_N_ROWS = 50  # markers kept per cluster, method and comparator role
W = Path("/exports/para-lipg-hpc/mdmanurung/cc-v8-work")
CAND = W / "tusoskill/cov95-2026-11-11/candidates"
RUN = W / "runs/v8-dev-aida-split-a-003"
PRED = {
    "celltypist": W / "mvp-aida-2026-11-10/predictions_celltypist.tsv",
    "popv": W / "mvp-aida-2026-11-10/predictions_popv.tsv",
    "singler": CAND / "c012_harmonizer/preds/aida_predictions_singler.tsv",
}
OUT = Path(__file__).resolve().parent / "demo"


def main():
    OUT.mkdir(exist_ok=True)

    # Cells: a stratified draw per cluster -------------------------------------------------
    adata = ad.read_h5ad(W / "prep-aida_split_a-2026-11-09/input.h5ad", backed="r")
    obs = adata.obs[["bench_cluster"]].astype({"bench_cluster": str})
    sizes = obs["bench_cluster"].value_counts()
    want = (sizes * TARGET_CELLS / len(obs)).round().astype(int).clip(lower=MIN_PER_CLUSTER)
    want = np.minimum(want, sizes)
    rng = np.random.default_rng(SEED)
    keep = []
    for cluster in sorted(sizes.index, key=lambda c: int(c)):
        members = np.flatnonzero(obs["bench_cluster"].to_numpy() == cluster)
        keep.append(rng.choice(members, size=int(want[cluster]), replace=False))
    keep = np.sort(np.concatenate(keep))  # original cell order
    barcodes = obs.index[keep].rename("barcode")
    small = ad.AnnData(
        X=np.zeros((len(keep), 0), dtype=np.float32),
        obs=obs.iloc[keep].copy(),
        obsm={"X_pca": np.asarray(adata.obsm["X_pca"])[keep].astype(np.float32)},
    )
    small.write_h5ad(OUT / "cells.h5ad")

    # Votes, tool predictions, truth: the same barcodes, by name ---------------------------
    votes = pd.read_csv(CAND / "c004_celllevel/cell_votes.tsv", sep="\t", dtype=str)
    votes = votes.set_index("barcode").loc[barcodes].reset_index()
    votes.to_csv(OUT / "cell_votes.tsv", sep="\t", index=False, lineterminator="\n")
    for tool, path in PRED.items():
        frame = pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str).set_index("barcode")
        frame.loc[barcodes].reset_index().to_csv(
            OUT / f"predictions_{tool}.tsv", sep="\t", index=False, lineterminator="\n"
        )
    truth = pd.read_csv(W / "prep-aida_split_a-2026-11-09/scoring/truth.tsv", sep="\t", dtype=str)
    truth = truth.set_index("barcode").loc[barcodes].reset_index()
    truth.to_csv(OUT / "truth.tsv", sep="\t", index=False, lineterminator="\n")

    # Vocabulary and marker table ------------------------------------------------------------
    vocab = pd.read_csv(RUN / "plan/vocabulary.tsv", sep="\t", dtype=str, keep_default_na=False)
    vocab.to_csv(OUT / "vocabulary.tsv", sep="\t", index=False, lineterminator="\n")
    markers = pd.read_csv(RUN / "annotation/levels/L1/markers.tsv", sep="\t", dtype=str, keep_default_na=False)
    markers = markers[markers["comparator_role"].eq("parent_remainder") & markers["rank"].astype(int).le(TOP_N_ROWS)]
    markers.to_csv(OUT / "markers.tsv", sep="\t", index=False, lineterminator="\n")
    print(f"demo: {len(keep)} cells, {obs.iloc[keep]['bench_cluster'].nunique()} clusters")


if __name__ == "__main__":
    main()
