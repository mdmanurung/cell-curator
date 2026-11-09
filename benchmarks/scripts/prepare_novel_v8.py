# Purpose: turn one protocol v8 novel dataset into the clustered input the annotators see.
#          Truth labels are written to a separate file that no annotator reads.
# Run:     sbatch benchmarks/scripts/prepare_novel_v8.sbatch DATASET RUN_DIR
#          DATASET is ts, indonesia or cross_tissue; RUN_DIR must not exist
# Inputs:  the dataset h5ad from acquire_v8.py (raw counts in raw.X)
# Outputs: RUN_DIR/input.h5ad     all genes, log1p CP10k in X, X_pca, obs: donor, tissue,
#                                 assay, bench_cluster; no truth label
#          RUN_DIR/scoring/truth.tsv  barcode, truth label, truth CL id (scoring only)
#          RUN_DIR/cluster_sizes.tsv  cells and donors per cluster
#          RUN_DIR/versions.txt

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# Parameters -----------------------------------------------------------------------
DOWNLOADS = Path("/exports/para-lipg-hpc/mdmanurung/cc-v8-work/acquire-2026-11-09/downloads")
FILES = {  # dataset_id from acquire_v8's downloads.tsv
    "ts": "c5d88abe-f23a-45fa-a534-788985e93dad.h5ad",
    "indonesia": "ca7d95ac-53aa-4a94-8c76-300c608ce6a0.h5ad",
    "cross_tissue": "1b9d8702-5af8-4142-85ed-020eb06ec4f6.h5ad",
    "aifi_dev": "e522d2cd-7927-4e59-a4ed-064009569279.h5ad",  # dev donors only
}
DEV_SPLIT = (
    Path(__file__).resolve().parents[2] / "benchmarks/data/aifi-hiha-dev-split/donor_split.tsv"
)
TS_TISSUES = ["blood", "spleen", "lymph node", "inguinal lymph node", "mesenteric lymph node"]
TS_ASSAYS = ["10x 3' v3", "10x 5' v2"]  # protocol v8, Q3
N_HVG = 2000  # Q17; seurat_v3 on raw counts, amendment 3
SCALE_CLIP = 10  # amendment 3
N_PCS = 50  # Q17
N_NEIGHBORS = 15  # Q17
RESOLUTION = 1.0  # Q15
SEED = 0  # Q17


def main():
    dataset, run_dir = sys.argv[1], Path(sys.argv[2])
    run_dir.mkdir(parents=True, exist_ok=False)  # a run dir is never overwritten
    (run_dir / "scoring").mkdir()

    # Raw counts, restricted as the protocol declares -----------------------------------
    if dataset == "aifi_dev":
        # backed read, then only the dev donors' rows; corpus donors are never loaded
        atlas = ad.read_h5ad(DOWNLOADS / FILES[dataset], backed="r")
        split = pd.read_csv(DEV_SPLIT, sep="\t")
        dev = set(split.loc[split["split"] == "dev", "donor_id"])
        rows = np.flatnonzero(atlas.obs["donor_id"].astype(str).isin(dev).to_numpy())
        data = ad.AnnData(
            X=atlas.raw.X[rows], obs=atlas.obs.iloc[rows].copy(), var=atlas.raw.var.copy()
        )
    else:
        atlas = ad.read_h5ad(DOWNLOADS / FILES[dataset])
        data = ad.AnnData(X=atlas.raw.X.copy(), obs=atlas.obs.copy(), var=atlas.raw.var.copy())
    # raw.var may lack symbols; take them from var by gene ID (fails as NaN if absent)
    data.var["feature_name"] = atlas.var["feature_name"].reindex(data.var_names).to_numpy()
    if data.var["feature_name"].isna().any():
        raise ValueError(f"{data.var['feature_name'].isna().sum()} genes have no symbol")
    if dataset == "ts":
        data = data[
            data.obs["tissue"].isin(TS_TISSUES) & data.obs["assay"].isin(TS_ASSAYS)
        ].copy()
    del atlas

    # One symbol per gene: the gene with the most counts keeps a shared symbol -------------
    symbol = data.var["feature_name"].astype(str)
    counts = pd.Series(np.asarray(data.X.sum(0)).ravel(), index=data.var_names)
    rank = counts.rank(ascending=False, method="first").to_numpy()
    first = (
        pd.Series(rank, index=data.var_names).groupby(symbol.to_numpy()).transform("min")
        == rank
    )
    data.var["gene"] = np.where(first, symbol, symbol + "--" + data.var_names)
    data.var_names = data.var["gene"].to_numpy()

    # HVG on raw counts, then log1p CP10k over all genes, scale the HVGs, PCA -----------
    sc.pp.highly_variable_genes(data, flavor="seurat_v3", n_top_genes=N_HVG)
    data.layers["counts"] = data.X.copy()  # the package's detection reads counts, as v6
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    hvg = data[:, data.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=SCALE_CLIP)
    sc.tl.pca(hvg, n_comps=N_PCS, random_state=SEED)
    sc.pp.neighbors(hvg, n_neighbors=N_NEIGHBORS, random_state=SEED)
    sc.tl.leiden(
        hvg, resolution=RESOLUTION, random_state=SEED, flavor="leidenalg", directed=False
    )
    data.obs["bench_cluster"] = hvg.obs["leiden"].astype(str).to_numpy()
    data.obsm["X_pca"] = hvg.obsm["X_pca"]

    # Truth goes to its own file; the annotator input carries none -------------------------
    data.obs[["cell_type", "cell_type_ontology_term_id"]].rename_axis(
        "barcode"
    ).reset_index().to_csv(run_dir / "scoring" / "truth.tsv", sep="\t", index=False)
    data.obs["annotation_scope"] = "ROOT"  # one scope; the package config needs the key
    keep = ["donor_id", "sex", "tissue", "assay", "bench_cluster", "annotation_scope"]
    out = ad.AnnData(
        X=data.X,
        obs=data.obs[keep].copy(),
        var=data.var[["gene"]],
        obsm=data.obsm,
        layers={"counts": data.layers["counts"]},
    )
    out.write_h5ad(run_dir / "input.h5ad")
    (
        out.obs.groupby("bench_cluster")
        .agg(cells=("donor_id", "size"), donors=("donor_id", "nunique"))
        .reset_index()
        .to_csv(run_dir / "cluster_sizes.tsv", sep="\t", index=False)
    )
    (run_dir / "versions.txt").write_text(
        f"scanpy {sc.__version__}\nanndata {ad.__version__}\n"
    )


if __name__ == "__main__":
    main()
