"""Cut the comparator training set out of the Human Immune Health Atlas.

Purpose: give CellTypist, SingleR, scANVI and popV the same labelled cells cell-curator
         learns from (protocol v8, Q6): the cells in sampled_cells.tsv, all genes.
Run:     python benchmarks/scripts/v8_prepare_training.py SAMPLED_TSV TRAIN_DIR [--check]
         (sbatch benchmarks/scripts/v8_prepare_training.sbatch)
Inputs:
    ATLAS_H5AD   raw counts in raw.X, obs cell_type and cell_type_ontology_term_id
    SAMPLED_TSV  learn_corpus_v8.py output: barcode, donor_id, cell_type,
                 cell_type_ontology_term_id (100 cells per donor and label, seed 1)
Outputs (TRAIN_DIR, must not exist):
    train.h5ad     X = log1p CP10k over all genes, layers["counts"], obs donor_id,
                   cell_type, cell_type_ontology_term_id (after the RELABEL swaps);
                   var_names = unique gene symbols
    label_map.tsv  label, cl_id: every training label with its Cell Ontology id
"""

from pathlib import Path

from v8_common import parse_args

# Parameters -----------------------------------------------------------------------
ATLAS_H5AD = Path(
    "/exports/para-lipg-hpc/mdmanurung/cc-v8-work/acquire-2026-11-09/downloads/"
    "e522d2cd-7927-4e59-a4ed-064009569279.h5ad"
)
TARGET_SUM = 1e4  # log1p CP10k, as the query inputs and v6 (CellTypist needs it too)
# Amendment 6, item 2 (user-approved): the two atlas state terms get the names and ids
# cell-curator's records use, so every method learns the same vocabulary.
RELABEL = {  # atlas cell_type: (training label, Cell Ontology id)
    "T cell": ("cycling T cell", "CL:4033069"),
    "CD8-positive, alpha-beta memory T cell": ("ISG+ memory CD8 T cell", "CL:0000909"),
}


def main():
    args = parse_args("training set for the v8 comparators", ["sampled_tsv"], ["train_dir"])
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc

    # Rows of the sampled cells, in file order ---------------------------------------------
    sampled = pd.read_csv(args.sampled_tsv, sep="\t", dtype=str)
    atlas = ad.read_h5ad(ATLAS_H5AD, backed="r")
    rows = atlas.obs_names.get_indexer(sampled["barcode"])
    if (rows < 0).any():
        raise ValueError(f"{(rows < 0).sum()} sampled barcodes are not in the atlas")
    order = np.argsort(rows)
    rows, sampled = rows[order], sampled.iloc[order].reset_index(drop=True)
    obs = atlas.obs.iloc[rows][["donor_id", "cell_type", "cell_type_ontology_term_id"]]
    obs = obs.astype(str)
    if not (obs["cell_type"].to_numpy() == sampled["cell_type"].to_numpy()).all():
        raise ValueError("atlas labels differ from sampled_cells.tsv")
    if not set(RELABEL) <= set(obs["cell_type"]):
        raise ValueError("a relabelled atlas term is absent from the sampled cells")
    obs = obs.copy()
    for old, (label, cl_id) in RELABEL.items():
        hit = obs["cell_type"] == old
        obs.loc[hit, ["cell_type", "cell_type_ontology_term_id"]] = [label, cl_id]

    # One symbol per gene: the gene with the most counts keeps a shared symbol ---------------
    counts = atlas.raw.X[rows]
    symbol = atlas.var["feature_name"].reindex(atlas.raw.var_names).astype(str)
    totals = pd.Series(np.asarray(counts.sum(0)).ravel(), index=atlas.raw.var_names)
    rank = totals.rank(ascending=False, method="first")
    first = rank.groupby(symbol.to_numpy()).transform("min") == rank
    names = np.where(first, symbol, symbol + "--" + atlas.raw.var_names)

    # Counts layer, then log1p CP10k in X -----------------------------------------------------
    data = ad.AnnData(
        X=counts.copy(), obs=obs, var=pd.DataFrame(index=pd.Index(names, name="gene"))
    )
    data.layers["counts"] = counts
    sc.pp.normalize_total(data, target_sum=TARGET_SUM)
    sc.pp.log1p(data)

    # Write ------------------------------------------------------------------------------------
    args.train_dir.mkdir(parents=True)
    data.write_h5ad(args.train_dir / "train.h5ad")
    labels = obs[["cell_type", "cell_type_ontology_term_id"]].drop_duplicates()
    labels.columns = ["label", "cl_id"]
    if not labels["label"].is_unique:
        raise ValueError("a training label carries more than one Cell Ontology id")
    labels.to_csv(args.train_dir / "label_map.tsv", sep="\t", index=False)
    print(f"{data.n_obs} cells, {data.n_vars} genes, {len(labels)} labels")


if __name__ == "__main__":
    main()
