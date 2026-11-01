"""Rank marker genes per benchmark cluster, the one input every LLM annotator receives.

Purpose: write, before any LLM annotator runs, the per-cluster marker table that
         GPTCelltype, mLLMCelltype, and CASSIA are all given, so they differ only in
         method and not in input.
Run:     uv run python benchmarks/data/aida-pbmc-protocol-v5/cluster_markers.py
Inputs:
    ../aida-phase1-v2-pbmc-hvg2000/aida_pbmc_clustered.h5ad   layers lognorm, obs bench_cluster
Outputs:
    cluster_markers.tsv   cluster, rank, gene, logfoldchange, pvals_adj (top 10 per cluster)
"""

from pathlib import Path

import anndata as ad
import scanpy as sc

HERE = Path(__file__).resolve().parent
INPUT = HERE.parent / "aida-phase1-v2-pbmc-hvg2000" / "aida_pbmc_clustered.h5ad"
OUT = HERE / "cluster_markers.tsv"
N_TOP = 10  # GPTCelltype's published default

if OUT.exists():
    raise ValueError(f"{OUT} exists; the v5 marker table is frozen")

adata = ad.read_h5ad(INPUT)
adata.X = adata.layers["lognorm"]
adata.var_names = adata.var["feature_name"].astype(str).to_numpy()
sc.tl.rank_genes_groups(adata, groupby="bench_cluster", method="wilcoxon", use_raw=False)
table = sc.get.rank_genes_groups_df(adata, group=None)
table = table.groupby("group", observed=True).head(N_TOP)
table["rank"] = table.groupby("group", observed=True).cumcount() + 1
table = table.rename(columns={"group": "cluster", "names": "gene", "logfoldchanges": "logfoldchange"})
table[["cluster", "rank", "gene", "logfoldchange", "pvals_adj"]].to_csv(OUT, sep="\t", index=False)
print(f"{table['cluster'].nunique()} clusters, {len(table)} rows")
