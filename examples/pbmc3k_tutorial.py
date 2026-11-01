"""Annotate the 10x PBMC 3k dataset with the shipped human PBMC marker panel.

Purpose: show cell-curator on real data, from raw counts to cluster labels, and put
         each label beside the cluster name from the scanpy tutorial.
Run:     uv run python examples/pbmc3k_tutorial.py [WORKDIR]
Inputs:
    sc.datasets.pbmc3k()            raw counts, 2,700 cells (downloaded, 5.9 MB)
    sc.datasets.pbmc3k_processed()  scanpy tutorial clusters, 2,638 cells (downloaded)
    skills/cell-curator/assets/markers/human-pbmc.json   shipped marker panel
    Cell Ontology 2026-06-08        fetched once by cell-curator, sha256-pinned
Outputs (under WORKDIR, a fresh temporary directory by default):
    results/pbmc3k-tutorial-v1/     the full run: evidence, labels, review packet
    printed table                   one row per cluster: scanpy name, label, confidence
"""

import logging
import sys
import tempfile
from pathlib import Path

import pandas as pd
import scanpy as sc

from cell_curator import CellCurator

ROOT = Path(__file__).resolve().parents[1]
WORKDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="pbmc3k-"))
# An installed package keeps the panel under <prefix>/share/cell-curator/assets/markers/.
MARKERS = ROOT / "skills" / "cell-curator" / "assets" / "markers" / "human-pbmc.json"
N_PCS = 50  # the scanpy tutorial's PCA size
SEED = 0

logging.getLogger("cell_curator").setLevel(logging.INFO)
logging.getLogger("cell_curator").addHandler(logging.StreamHandler())
(WORKDIR / "data").mkdir(parents=True, exist_ok=True)
sc.settings.datasetdir = WORKDIR / "data"

# Load counts and clusters --------------------------------------------------------
# Counts come from the raw dataset; clusters come from the scanpy tutorial, which
# filtered cells, so the join keeps only the cells it kept.

adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
clusters = sc.datasets.pbmc3k_processed().obs["louvain"]
n_raw = adata.n_obs
adata = adata[adata.obs_names.isin(clusters.index)].copy()
print(f"kept {adata.n_obs} of {n_raw} cells that the scanpy tutorial clustered")

names = clusters.loc[adata.obs_names]
if not names.index.equals(adata.obs_names):
    raise ValueError("cluster labels are not aligned to cells")
# Cluster IDs are codes, so cell-curator never sees the tutorial's cell type names.
adata.obs["cluster"] = names.cat.codes.astype(str).to_numpy()
adata.obs["scanpy_name"] = names.to_numpy()

# Normalize and embed -------------------------------------------------------------

# The scanpy tutorial's gene filter. It also cuts the CPU marker-ranking time, which
# scales with the number of genes.
n_genes = adata.n_vars
sc.pp.filter_genes(adata, min_cells=3)
print(f"kept {adata.n_vars} of {n_genes} genes detected in at least 3 cells")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
adata.var["gene_symbol"] = adata.var_names
sc.pp.pca(adata, n_comps=N_PCS, random_state=SEED)

# Annotate ------------------------------------------------------------------------
# preauthorized=True declares the Gate A approval up front for a scripted run; an
# interactive run stops and asks a person instead.

curator = CellCurator(
    adata=adata,
    organism="human",
    tissue="peripheral blood",
    cluster_key="cluster",
    markers=MARKERS,
    run_id="pbmc3k-tutorial-v1",
    output_root=WORKDIR / "results",
    mode="annotation",
    embedding_key="X_pca",
    developmental_stage="adult",
    condition="healthy",
    experimental_context="10x Genomics PBMC 3k, scanpy tutorial clustering",
    preauthorized=True,
    reviewer="pbmc3k-tutorial",
    assumptions=["The shipped 13-type human PBMC panel defines the L1 vocabulary."],
    overrides={
        # Overrides replace lists, so the marker provider is listed again beside the
        # ontology. online only so the ontology can be fetched once; nothing else is sent.
        "knowledge": {
            "network": "online",
            "providers": [
                {"kind": "local_markers", "path": str(MARKERS), "version": "human-pbmc-v1"},
                {
                    "kind": "remote_ontology",
                    "url": "http://purl.obolibrary.org/obo/cl/releases/2026-06-08/cl.obo",
                    "cache_path": str(WORKDIR / "cache" / "cl-2026-06-08.obo"),
                    "sha256": "7e19b4aef8e7fe7720b0c4474c6a5b8713c904ea9ab7b7d5d96f2cb28b8fc51f",
                    "version": "2026-06-08",
                    "timeout_seconds": 120,
                },
            ],
        }
    },
)
curator.inspect()
curator.prepare()
curator.audit_parents()
curator.evidence()
curator.annotate()
curator.validate_run()

# Compare with the scanpy tutorial ------------------------------------------------

labels = pd.read_csv(
    WORKDIR / "results" / "pbmc3k-tutorial-v1" / "labels" / "L1.tsv",
    sep="\t",
    dtype={"cluster_id": str},
)
names_by_cluster = (
    adata.obs[["cluster", "scanpy_name"]]
    .drop_duplicates()
    .rename(columns={"cluster": "cluster_id"})
)
table = labels.merge(names_by_cluster, on="cluster_id", how="left", validate="1:1")
print(table[["cluster_id", "scanpy_name", "label", "cl_id", "confidence"]].to_string(index=False))
print(f"run directory: {WORKDIR / 'results' / 'pbmc3k-tutorial-v1'}")
