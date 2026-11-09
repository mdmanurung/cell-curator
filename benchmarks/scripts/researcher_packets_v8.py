"""Write the evidence packet the three v8 researcher agents read for one novel dataset.

Purpose: protocol v8 comparator and consistency measure. Each researcher annotates every
cluster of a novel dataset from the same packet. The packet holds measured numbers only:
no truth label, no cell-curator call, no dataset-specific hint beyond the obs columns.
Run:     .venv/bin/python benchmarks/scripts/researcher_packets_v8.py DATASET INPUT_H5AD OUT_DIR
         DATASET is a label such as ts, indonesia or cross_tissue; OUT_DIR must not exist
Inputs:  INPUT_H5AD  prepare_novel_v8.py output: log1p CP10k in X, gene symbols as var_names,
                     obs donor, tissue, assay, bench_cluster [prepare_novel_v8.sbatch]
Outputs: OUT_DIR/packet.json  the packet, one entry per cluster (the researchers' input)
         OUT_DIR/packet.md    the same numbers as tables, for reading
         OUT_DIR/versions.txt scanpy and package versions
"""

import json
import sys
from importlib import metadata
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# Parameters -----------------------------------------------------------------------
N_TOP_GENES = 20  # user decision 2026-11-09: top genes per cluster by Wilcoxon vs rest
RANK_METHOD = "wilcoxon"  # cell-level, cluster versus all other cells; same as the user asked
TOP_DONORS_SHOWN = 3  # donors listed per cluster with their share of the cluster
# Canonical lineage panel, proposed (see "Open questions" in benchmarks/v8-researcher/mapping_rule.md). Fixed before any packet is written; the same list goes to all datasets.
# A gene missing from var_names is listed in the packet as not measured, not dropped.
LINEAGE_PANEL = {
    "pan-leukocyte": ["PTPRC"],
    "T cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "IL7R", "CCR7", "SELL", "TCF7"],
    "T effector, memory, regulatory": [
        "GZMK",
        "GZMB",
        "GZMH",
        "PRF1",
        "CCL5",
        "FOXP3",
        "IL2RA",
        "CTLA4",
        "KLRB1",
        "TRDC",
        "TRGC1",
    ],  # fmt: skip
    "NK and ILC": ["NKG7", "GNLY", "KLRD1", "NCAM1", "FCGR3A", "KLRC1"],
    "B cell": ["CD79A", "MS4A1", "CD19", "CD27", "IGHM", "IGHD", "TCL1A", "CXCR5"],
    "plasma cell": ["JCHAIN", "MZB1", "SDC1", "IGHG1", "IGHA1"],
    "monocyte and macrophage": [
        "CD14",
        "LYZ",
        "S100A8",
        "FCN1",
        "CD68",
        "C1QA",
        "MARCO",
        "MRC1",
    ],  # fmt: skip
    "dendritic cell": ["CST3", "CLEC9A", "CD1C", "FCER1A", "LILRA4", "CLEC4C", "IL3RA"],
    "granulocyte and mast": ["CSF3R", "FCGR3B", "CPA3", "KIT"],
    "platelet, erythroid, progenitor": ["PPBP", "PF4", "HBB", "HBA1", "CD34"],
    "proliferation": ["MKI67", "TOP2A"],
    "tissue residence and activation": ["CD69", "ITGAE", "PDCD1", "HLA-DRA"],
    "non-immune": ["EPCAM", "PECAM1", "COL1A1"],
}
PANEL_GENES = [gene for genes in LINEAGE_PANEL.values() for gene in genes]


# Read and check ------------------------------------------------------------------------


def main():
    dataset, h5ad, out_dir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    out_dir.mkdir(parents=True)  # fails if the run dir exists

    adata = ad.read_h5ad(h5ad)
    obs = adata.obs
    # Truth leaks through any extra obs column; the packet reads only these.
    obs = obs[["donor_id", "tissue", "assay", "bench_cluster"]].copy()
    obs["bench_cluster"] = obs["bench_cluster"].astype(str)
    adata.obs = obs
    adata.obs["bench_cluster"] = adata.obs["bench_cluster"].astype("category")
    if not adata.var_names.is_unique:
        raise ValueError("var_names (gene symbols) must be unique")

    # Wilcoxon, each cluster versus the rest ---------------------------------------------

    # scanpy ranks on X as given; X is log1p CP10k so no further normalisation is applied.
    sc.tl.rank_genes_groups(
        adata, "bench_cluster", method=RANK_METHOD, reference="rest", use_raw=False
    )
    ranked = sc.get.rank_genes_groups_df(adata, group=None)
    ranked = ranked.sort_values(["group", "scores"], ascending=[True, False])
    ranked = ranked.groupby("group", observed=True).head(N_TOP_GENES)

    # Detection and mean expression, inside and outside each cluster --------------------

    clusters = list(adata.obs["bench_cluster"].cat.categories)
    top_genes = sorted(set(ranked["names"]) | set(PANEL_GENES) & set(adata.var_names))
    expr = adata[:, top_genes].X
    expr = expr.toarray() if hasattr(expr, "toarray") else np.asarray(expr)
    expr = pd.DataFrame(expr, index=adata.obs_names, columns=top_genes)
    in_cluster = {c: (adata.obs["bench_cluster"] == c).to_numpy() for c in clusters}
    detect_in = pd.DataFrame({c: (expr[in_cluster[c]] > 0).mean() for c in clusters})
    detect_out = pd.DataFrame({c: (expr[~in_cluster[c]] > 0).mean() for c in clusters})
    mean_in = pd.DataFrame({c: expr[in_cluster[c]].mean() for c in clusters})

    # One entry per cluster --------------------------------------------------------------

    def mix(cluster, column):
        share = adata.obs.loc[in_cluster[cluster], column].value_counts(normalize=True)
        return {str(k): round(float(v), 4) for k, v in share[share > 0].items()}

    packet = {
        "dataset": dataset,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "expression_scale": "log1p of counts per 10,000",
        "ranking": f"{RANK_METHOD}, cell level, each cluster versus all other cells, "
        f"top {N_TOP_GENES} by test score, positive direction",
        "panel_genes_not_measured": [g for g in PANEL_GENES if g not in adata.var_names],
        "clusters": [],
    }
    for cluster in clusters:
        donors = adata.obs.loc[in_cluster[cluster], "donor_id"].value_counts(normalize=True)
        top = ranked[ranked["group"] == cluster]
        packet["clusters"].append(
            {
                "cluster": cluster,
                "n_cells": int(in_cluster[cluster].sum()),
                "n_donors": int((donors > 0).sum()),
                "top_donors": {
                    str(k): round(float(v), 4) for k, v in donors.head(TOP_DONORS_SHOWN).items()
                },
                "tissue_mix": mix(cluster, "tissue"),
                "assay_mix": mix(cluster, "assay"),
                "top_genes": [
                    {
                        "gene": row.names,
                        "log2_fold_change": round(float(row.logfoldchanges), 3),
                        "detected_inside": round(float(detect_in.loc[row.names, cluster]), 4),
                        "detected_outside": round(float(detect_out.loc[row.names, cluster]), 4),
                        "adjusted_p": float(f"{row.pvals_adj:.3g}"),
                    }
                    for row in top.itertuples()
                ],
                "lineage_panel": {
                    gene: {
                        "mean_log1p_cp10k": round(float(mean_in.loc[gene, cluster]), 3),
                        "detected_inside": round(float(detect_in.loc[gene, cluster]), 4),
                    }
                    for gene in PANEL_GENES
                    if gene in adata.var_names
                },
            }
        )

    # Write ------------------------------------------------------------------------------

    (out_dir / "packet.json").write_text(json.dumps(packet, indent=1))
    (out_dir / "packet.md").write_text(render_markdown(packet))
    (out_dir / "versions.txt").write_text(
        "\n".join(f"{n} {metadata.version(n)}" for n in ("scanpy", "anndata", "pandas")) + "\n"
    )


def render_markdown(packet):
    lines = [
        f"# Cluster packet: {packet['dataset']}",
        "",
        f"{packet['n_cells']} cells, {packet['n_genes']} genes. Expression is "
        f"{packet['expression_scale']}. Genes ranked by {packet['ranking']}.",
    ]
    if packet["panel_genes_not_measured"]:
        lines += [
            "",
            "Panel genes not measured: " + ", ".join(packet["panel_genes_not_measured"]),
        ]
    for c in packet["clusters"]:
        lines += [
            "",
            f"## Cluster {c['cluster']}",
            "",
            f"{c['n_cells']} cells, {c['n_donors']} donors. Largest donor shares: {c['top_donors']}.",
            f"Tissue mix: {c['tissue_mix']}. Assay mix: {c['assay_mix']}.",
            "",
            "| gene | log2 FC | detected inside | detected outside | adjusted p |",
            "| --- | --- | --- | --- | --- |",
        ]
        lines += [
            f"| {g['gene']} | {g['log2_fold_change']} | {g['detected_inside']} "
            f"| {g['detected_outside']} | {g['adjusted_p']} |"
            for g in c["top_genes"]
        ]
        lines += [
            "",
            "| panel gene | mean expression | detected inside |",
            "| --- | --- | --- |",
        ]
        lines += [
            f"| {gene} | {v['mean_log1p_cp10k']} | {v['detected_inside']} |"
            for gene, v in c["lineage_panel"].items()
        ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
