"""Write a benchmark slice in the form CellTypist reads: log1p CP10k in X, symbols as names.

Purpose: CellTypist annotates `adata.X` by gene symbol; the benchmark inputs keep the
         log-normalized values in `layers["lognorm"]` and, for v6, Ensembl IDs as names.
Run:     python benchmarks/scripts/prepare_celltypist_input.py INPUT_H5AD OUT_H5AD [SYMBOLS_TSV]
Inputs:
    INPUT_H5AD    v5 aida_pbmc_clustered.h5ad or v6 aida_pbmc_wholegenome_clustered.h5ad
    SYMBOLS_TSV   v6 gene_symbols.tsv; omitted, var feature_name is used (v5)
Outputs:
    OUT_H5AD      X = lognorm, var_names = gene symbols, obs bench_cluster kept
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd

INPUT, OUT = Path(sys.argv[1]), Path(sys.argv[2])
SYMBOLS = Path(sys.argv[3]) if len(sys.argv) > 3 else None
if OUT.exists():
    raise ValueError(f"refusing to overwrite {OUT}")

adata = ad.read_h5ad(INPUT)
if SYMBOLS is None:
    names = adata.var["feature_name"].astype(str).to_numpy()
else:
    names = pd.read_csv(SYMBOLS, sep="\t").set_index("feature_id")["gene_symbol"]
    names = names.loc[adata.var_names].to_numpy()
out = ad.AnnData(
    X=adata.layers["lognorm"].copy(),
    obs=adata.obs[["bench_cluster"]].copy(),
    var=pd.DataFrame(index=pd.Index(names, name="gene_symbol")),
)
out.write_h5ad(OUT)
print(f"{out.n_obs} cells x {out.n_vars} genes")
