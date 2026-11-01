"""Split the AIDA slice's donors in half for the reference-trained comparators.

Purpose: fix, before any v5 comparator runs, which donors train SingleR, scANVI, and
         popV (split A) and which donors every system is scored on (split B).
Run:     uv run python benchmarks/data/aida-pbmc-protocol-v5/split_donors.py
Inputs:
    ../aida-phase1-v2-pbmc-hvg2000/aida_pbmc_clustered.h5ad   obs donor_id (pinned v5 input)
Outputs:
    donor_split.tsv    donor_id, split (A or B), n_cells
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
INPUT = HERE.parent / "aida-phase1-v2-pbmc-hvg2000" / "aida_pbmc_clustered.h5ad"
OUT = HERE / "donor_split.tsv"
SEED = 0

if OUT.exists():
    raise ValueError(f"{OUT} exists; the v5 donor split is frozen")

donors = ad.read_h5ad(INPUT, backed="r").obs["donor_id"].astype(str)
counts = donors.value_counts().sort_index()
order = np.random.default_rng(SEED).permutation(len(counts))
split = np.where(np.arange(len(counts)) < len(counts) // 2, "A", "B")
table = pd.DataFrame(
    {"donor_id": counts.index[order], "split": split, "n_cells": counts.to_numpy()[order]}
).sort_values("donor_id")
table.to_csv(OUT, sep="\t", index=False)
print(table.groupby("split").agg(donors=("donor_id", "size"), cells=("n_cells", "sum")))
