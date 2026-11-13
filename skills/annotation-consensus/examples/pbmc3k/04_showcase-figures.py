# Purpose: draw the PBMC3K showcase figures from the two consensus runs and write the cluster table.
# Run:     .venv/bin/python 04_showcase-figures.py RUN_DIR      (RUN_DIR/figures must not exist)
# Inputs:  RUN_DIR/pbmc3k.h5ad                         [01_prep-pbmc3k.py]
#          RUN_DIR/out_h5/{consensus,calls}.tsv        run 1: rule h5, marker co-voter [annotation_consensus.py]
#          RUN_DIR/out_majority/{consensus,calls}.tsv  run 2: plain majority, no marker vote
# Outputs: RUN_DIR/figures/fig1_umap.png          UMAP: scanpy name, strict label, best-guess label
#          RUN_DIR/figures/fig2_tool_votes.png    per cluster: each tool's call and its share of the cluster
#          RUN_DIR/figures/fig3_runs.png          per cluster: cells called by run 1 and run 2
#          RUN_DIR/figures/fig4_crosstab.png      consensus label (run 1 best guess) against scanpy name
#          RUN_DIR/figures/fig5_markers.png       vocabulary markers of each best-guess label, per cluster
#          RUN_DIR/figures/cluster_table.tsv      one row per cluster: scanpy name, both runs' labels
# PBMC3K is one donor, so every point is a cell; the scanpy name is a reference, not truth.

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
run_dir = Path(sys.argv[1])
fig_dir = run_dir / "figures"
fig_dir.mkdir()

# Parameters -----------------------------------------------------------------------
SEED = 0
N_MARKERS = 3  # positive markers shown per best-guess label
MM = 1 / 25.4
plt.rcParams.update({"font.size": 7, "axes.linewidth": 0.5, "figure.dpi": 150})

# Load ------------------------------------------------------------------------------
adata = ad.read_h5ad(run_dir / "pbmc3k.h5ad")
read = lambda name, file: pd.read_csv(  # noqa: E731
    run_dir / name / file, sep="\t", keep_default_na=False, dtype=str
)
cons = {
    n: read(f"out_{n}", "consensus.tsv").set_index("cluster_id") for n in ("h5", "majority")
}
calls = read("out_h5", "calls.tsv").set_index("barcode").loc[adata.obs_names]
vocab = pd.read_csv(HERE / "vocabulary.tsv", sep="\t", keep_default_na=False)
if not calls["cluster_id"].equals(adata.obs["cluster"].astype(str)):
    raise ValueError("calls.tsv clusters differ from the h5ad clusters")

names = adata.obs.drop_duplicates("cluster").set_index("cluster")["scanpy_name"]
order = sorted(names.index, key=int)
short = lambda s: s.replace(", alpha-beta", "").replace("thymus-derived ", "")  # noqa: E731

# Cluster table ---------------------------------------------------------------------
table = pd.DataFrame(
    {
        "scanpy_name": names.loc[order],
        "run1_h5_strict": cons["h5"].loc[order, "label"],
        "run1_h5_best_guess": cons["h5"].loc[order, "label_best_guess"],
        "run2_majority_strict": cons["majority"].loc[order, "label"],
        "run2_majority_best_guess": cons["majority"].loc[order, "label_best_guess"],
        "n_cells": cons["h5"].loc[order, "n_cells"],
    }
)
table.to_csv(fig_dir / "cluster_table.tsv", sep="\t")

# Fig 1: UMAP ------------------------------------------------------------------------
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
adata.obs["strict"] = calls["label_strict"].map(short).to_numpy()
adata.obs["best_guess"] = calls["label_best_guess"].map(short).to_numpy()
fig, axes = plt.subplots(3, 1, figsize=(170 * MM, 230 * MM))
for ax, key, title in zip(
    axes,
    ("scanpy_name", "strict", "best_guess"),
    (
        "scanpy tutorial name (reference)",
        "consensus strict (run 1)",
        "consensus best guess (run 1)",
    ),
    strict=True,
):
    sc.pl.umap(
        adata,
        color=key,
        ax=ax,
        show=False,
        title=title,
        size=5,
        legend_fontsize=5,
        frameon=False,
    )
    ax.legend(
        loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=5, markerscale=0.6
    )
fig.tight_layout()
fig.savefig(fig_dir / "fig1_umap.png", dpi=200)
plt.close(fig)

# Fig 2: tool votes per cluster ------------------------------------------------------
tools = ("celltypist", "popv", "singler")
fig, axes = plt.subplots(1, 3, figsize=(260 * MM, 70 * MM), sharey=True)
for ax, tool in zip(axes, tools, strict=True):
    share = np.array(
        [json.loads(cons["h5"].loc[c, "tool_calls"])[tool]["cluster_share"] for c in order]
    )
    ax.barh(range(len(order)), share.astype(float), color="0.7")
    for i, c in enumerate(order):
        ax.text(
            0.01,
            i,
            short(json.loads(cons["h5"].loc[c, "tool_calls"])[tool]["label"]),
            va="center",
            fontsize=5.5,
        )
    ax.set_xlim(0, 1)
    ax.set_title(tool)
    ax.set_xlabel("share of the cluster's cells")
axes[0].set_yticks(range(len(order)), [f"{c} {names[c]}" for c in order])
axes[0].invert_yaxis()
fig.tight_layout()
fig.savefig(fig_dir / "fig2_tool_votes.png", dpi=200)
plt.close(fig)

# Fig 3: share of each cluster's cells with a strict label, per run ---------------------------------
share = {
    n: (
        read(f"out_{n}", "calls.tsv")
        .assign(called=lambda d: d["label_strict"].ne("Unknown"))
        .groupby("cluster_id")["called"]
        .mean()
        .loc[order]
    )
    for n in cons
}
fig, ax = plt.subplots(figsize=(130 * MM, 75 * MM))
width = 0.38
for k, (n, label) in enumerate(
    (("h5", "run 1: tools + marker vote (h5)"), ("majority", "run 2: tools only (majority)"))
):
    ax.bar(np.arange(len(order)) + (k - 0.5) * width, share[n].to_numpy(), width, label=label)
ax.set_xticks(range(len(order)), [names[c] for c in order], rotation=40, ha="right")
ax.set_ylim(0, 1.25)
ax.set_ylabel("cells with a strict label")
ax.legend(frameon=False, fontsize=6, loc="upper right")
fig.tight_layout()
fig.savefig(fig_dir / "fig3_runs.png", dpi=200)
plt.close(fig)

# Fig 4: consensus label against scanpy name ---------------------------------------------
cross = pd.crosstab(adata.obs["best_guess"], adata.obs["scanpy_name"], normalize="columns")
fig, ax = plt.subplots(figsize=(150 * MM, 90 * MM))
im = ax.imshow(cross.to_numpy(), aspect="auto", cmap="Greys", vmin=0, vmax=1)
ax.set_xticks(range(cross.shape[1]), cross.columns, rotation=40, ha="right")
ax.set_yticks(range(cross.shape[0]), cross.index)
fig.colorbar(im, ax=ax, label="share of the scanpy cluster")
fig.tight_layout()
fig.savefig(fig_dir / "fig4_crosstab.png", dpi=200)
plt.close(fig)

# Fig 5: vocabulary markers of each best-guess label, per cluster --------------------------
labels = pd.unique(cons["h5"].loc[order, "label_best_guess"])
genes = {}
for lab in labels:
    row = vocab[vocab["label"].eq(lab) & vocab["axis"].eq("cell_type")]
    genes[lab] = [
        g for g in row["positive_markers"].iloc[0].split("|") if g in adata.var_names
    ][:N_MARKERS]
var_names = {short(k): v for k, v in genes.items() if v}
adata.obs["cluster_name"] = pd.Categorical(
    adata.obs["cluster"].map(lambda c: f"{c} {names[c]}"),
    categories=[f"{c} {names[c]}" for c in order],
)
dp = sc.pl.dotplot(
    adata, var_names, groupby="cluster_name", layer="logcounts", show=False, return_fig=True
)
dp.savefig(fig_dir / "fig5_markers.png", dpi=200, bbox_inches="tight")
plt.close("all")
sc.logging.print_versions()
