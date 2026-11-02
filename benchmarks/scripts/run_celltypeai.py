"""Run CellTypeAI on the v5 benchmark clusters through a local Ollama server.

Purpose: the CellTypeAI comparator for protocol v5 (Daw et al., Bioinformatics 2026).
         It calls CellTypeAI's own ranking and prompt functions in the loop
         `cell_annotator` uses, so each repeat's answer is kept for the declared
         confidence (frequency of the modal answer), which `cell_annotator` discards.
Run:     python benchmarks/scripts/run_celltypeai.py OUT_DIR   (Ollama on localhost:11434)
Inputs:
    benchmarks/data/aida-phase1-v2-pbmc-hvg2000/aida_pbmc_clustered.h5ad  layers lognorm,
                                                                          obs bench_cluster
Outputs (OUT_DIR, must not exist):
    cell_predictions.tsv   barcode, predicted_label, uncertainty (1 - modal frequency)
    cluster_answers.tsv    cluster, repeat, answer
    cluster_genes.tsv      cluster, rank, gene (CellTypeAI's own Wilcoxon ranking)
"""

import sys
from collections import Counter
from pathlib import Path

import anndata as ad
import pandas as pd
from celltypeai.h5ad_to_json import _rank_var
from celltypeai.prompt_submit import _cell_type_prompt

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "benchmarks/data/aida-phase1-v2-pbmc-hvg2000/aida_pbmc_clustered.h5ad"
OUT = Path(sys.argv[1])
# CellTypeAI's defaults; "Blood" is in its tissue list.
SPECIES, TISSUE, MODEL, NUM_GENES, N_ITERATIONS = "Human", "Blood", "phi4:14b", 200, 3

OUT.mkdir(parents=True)
adata = ad.read_h5ad(INPUT)
adata.X = adata.layers["lognorm"]
adata.var_names = adata.var["feature_name"].astype(str).to_numpy()
adata.obs["louvain"] = adata.obs["bench_cluster"].astype(str).astype("category")

# Rank and prompt ----------------------------------------------------------------------
# Answers are keyed by group position; the position's cluster name comes from the
# ranking itself, so no answer is attached to the wrong cluster.

gene_data, num_groups = _rank_var(adata, "louvain")
clusters = list(gene_data.dtype.names)
answers = [
    _cell_type_prompt(
        SPECIES, TISSUE, adata, MODEL, gene_data, num_groups, NUM_GENES, False
    )[0]
    for _ in range(N_ITERATIONS)
]
long = pd.DataFrame(
    [
        {"cluster": clusters[i], "repeat": r, "answer": answer.get(str(i), "")}
        for r, answer in enumerate(answers)
        for i in range(num_groups)
    ]
)
long.to_csv(OUT / "cluster_answers.tsv", sep="\t", index=False)
pd.DataFrame(
    [
        {"cluster": name, "rank": rank + 1, "gene": gene_data[rank][name]}
        for name in clusters
        for rank in range(NUM_GENES)
    ]
).to_csv(OUT / "cluster_genes.tsv", sep="\t", index=False)

# Mode and confidence -------------------------------------------------------------------

modal = {}
for name, block in long.groupby("cluster"):
    label, count = Counter(block["answer"]).most_common(1)[0]
    modal[name] = (label or "Unknown", count / N_ITERATIONS)
cells = adata.obs["louvain"].astype(str)
pd.DataFrame(
    {
        "barcode": adata.obs_names,
        "predicted_label": cells.map(lambda c: modal[c][0]).to_numpy(),
        "uncertainty": cells.map(lambda c: 1 - modal[c][1]).to_numpy(),
    }
).to_csv(OUT / "cell_predictions.tsv", sep="\t", index=False)
print(pd.Series({c: v[0] for c, v in modal.items()}).to_string())
