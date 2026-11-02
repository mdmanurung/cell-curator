"""Write the v5 donor splits as h5ad and Matrix Market for the reference-trained comparators.

Purpose: give SingleR (R), scANVI, and popV the same split A (train) and split B (score)
         cells, from the frozen donor split.
Run:     uv run python benchmarks/scripts/export_splits.py OUT_DIR [INPUT_H5AD SYMBOLS_TSV]
         (protocol v6 passes the whole-transcriptome input and its gene_symbols.tsv)
Inputs:
    benchmarks/data/aida-phase1-v2-pbmc-hvg2000/aida_pbmc_clustered.h5ad   pinned v5 input
    benchmarks/data/aida-pbmc-protocol-v5/donor_split.tsv                   frozen split
Outputs (OUT_DIR, must not exist):
    A.h5ad, B.h5ad            the two splits, all layers and obs
    {A,B}/lognorm.mtx         genes x cells log1p CP10k, for R
    {A,B}/genes.tsv           gene symbols, one per line
    {A,B}/barcodes.tsv        cell barcodes, one per line
    A/labels.tsv              barcode, cell_type (AIDA truth, for training only)
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import io

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "benchmarks/data/aida-phase1-v2-pbmc-hvg2000/aida_pbmc_clustered.h5ad"
SPLIT = ROOT / "benchmarks/data/aida-pbmc-protocol-v5/donor_split.tsv"
OUT = Path(sys.argv[1])
if len(sys.argv) > 2:
    INPUT = Path(sys.argv[2])
SYMBOLS = Path(sys.argv[3]) if len(sys.argv) > 3 else None

OUT.mkdir(parents=True)
adata = ad.read_h5ad(INPUT)
split = pd.read_csv(SPLIT, sep="\t", dtype={"donor_id": str})
adata.obs = adata.obs.merge(
    split[["donor_id", "split"]], on="donor_id", how="left", validate="m:1"
).set_index(adata.obs.index)
if adata.obs["split"].isna().any():
    raise ValueError("cells with a donor missing from the split")
if SYMBOLS is None:
    adata.var_names = adata.var["feature_name"].astype(str).to_numpy()
else:
    symbols = pd.read_csv(SYMBOLS, sep="\t").set_index("feature_id")["gene_symbol"]
    adata.var_names = symbols.loc[adata.var_names].to_numpy()

for name in ("A", "B"):
    part = adata[adata.obs["split"].eq(name)].copy()
    part.write_h5ad(OUT / f"{name}.h5ad")
    folder = OUT / name
    folder.mkdir()
    io.mmwrite(folder / "lognorm.mtx", part.layers["lognorm"].T.tocoo())
    pd.Series(part.var_names).to_csv(folder / "genes.tsv", index=False, header=False)
    pd.Series(part.obs_names).to_csv(folder / "barcodes.tsv", index=False, header=False)
    if name == "A":
        part.obs[["cell_type"]].rename_axis("barcode").to_csv(folder / "labels.tsv", sep="\t")
    print(f"{name}: {part.n_obs} cells, {part.obs['donor_id'].nunique()} donors")
