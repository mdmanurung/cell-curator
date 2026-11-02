"""Give every gene in the whole-transcriptome slice one unique symbol, for protocol v6.

Purpose: the engine needs unique feature names, and 1,469 Census genes share a symbol.
         For each shared symbol, the gene with the most total counts keeps it; the
         others become SYMBOL--ENSEMBL_ID, so no marker can match a low-count copy.
Run:     python benchmarks/data/aida-phase1-v2-pbmc-wholegenome/gene_symbols.py
Inputs:
    aida_pbmc_wholegenome_clustered.h5ad   [acquire_wholegenome.py]; var feature_name,
                                           layers counts
Outputs:
    gene_symbols.tsv   feature_id, gene_symbol (one row per gene, symbols unique)
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "gene_symbols.tsv"
if OUT.exists():
    raise ValueError(f"{OUT} exists; the v6 symbol table is frozen")

adata = ad.read_h5ad(HERE / "aida_pbmc_wholegenome_clustered.h5ad")
genes = pd.DataFrame(
    {
        "feature_id": adata.var_names.astype(str),
        "symbol": adata.var["feature_name"].astype(str).to_numpy(),
        "total": np.asarray(adata.layers["counts"].sum(axis=0)).ravel(),
    }
)
ranked = genes.sort_values(["symbol", "total", "feature_id"], ascending=[True, False, True])
keeps = ~ranked.duplicated("symbol")
ranked["gene_symbol"] = np.where(keeps, ranked["symbol"], ranked["symbol"] + "--" + ranked["feature_id"])
table = genes[["feature_id"]].merge(
    ranked[["feature_id", "gene_symbol"]], on="feature_id", how="left", validate="1:1"
)
if table["gene_symbol"].duplicated().any():
    raise ValueError("symbols are still not unique")
table.to_csv(OUT, sep="\t", index=False)
print(f"{len(table)} genes; {int((~keeps).sum())} renamed with their Ensembl ID")
