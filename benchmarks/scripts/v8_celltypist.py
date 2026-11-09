"""Train CellTypist on the v8 training cells and annotate a v8 query.

Purpose: the CellTypist comparator for protocol v8 (Q6): a model retrained on the AIFI
         corpus cells, not the pretrained Immune_All_Low.pkl. Prediction follows the v6
         lane: per cell and majority vote over bench_cluster, uncertainty = 1 - conf_score
         of the per-cell label. CellTypist cannot abstain here: the default
         mode="best match" always assigns a label.
Run:     python benchmarks/scripts/v8_celltypist.py TRAIN_DIR QUERY_H5AD OUT_DIR [--check]
         (sbatch benchmarks/scripts/v8_celltypist.sbatch; env envs/cc-celltypist)
Inputs:
    TRAIN_DIR/train.h5ad   [v8_prepare_training.py]
    QUERY_H5AD             input.h5ad: log1p CP10k in X, symbols as var_names, obs bench_cluster
Outputs (OUT_DIR, must not exist):
    model.pkl                       the trained CellTypist model
    per_cell/cell_predictions.tsv   barcode, predicted_label, uncertainty
    majority/cell_predictions.tsv   same, predicted_label = majority-vote label
"""

import os

from v8_common import parse_args

# Parameters -----------------------------------------------------------------------
# OPEN QUESTION: no v6 value (v6 used the pretrained model). Provisional: the route the
# CellTypist docs give for large training sets (SGD, with a first pass that keeps the
# top genes per label); the package default (False, False) densifies 200k x 32k.
USE_SGD = True
FEATURE_SELECTION = True
TOP_GENES = 300  # package default, top genes per label kept after the first pass
BALANCE_CELL_TYPE = False  # package default; the sampling is already capped per donor and label
EPOCHS = 10  # package default
N_JOBS = 4  # speed only (8 on the 8-CPU sbatch; 4 in the 4-CPU allocation)
SEED = 0  # as the other comparators
# A model.pkl from an earlier run of this script on the same TRAIN_DIR; skips training.
# The fit does not see the query, so the model is the same one.
REUSE_MODEL = os.environ.get("REUSE_MODEL", "")
CLUSTER_COLUMN = "bench_cluster"  # v6: majority vote over the benchmark clusters


def main():
    args = parse_args(
        "CellTypist retrained on the v8 corpus cells", ["train_dir", "query_h5ad"], ["out_dir"]
    )
    import anndata as ad
    import celltypist

    query = ad.read_h5ad(args.query_h5ad)

    # Train -------------------------------------------------------------------------------------
    args.out_dir.mkdir(parents=True)
    if REUSE_MODEL:
        model = celltypist.models.Model.load(REUSE_MODEL)
    else:
        train = ad.read_h5ad(args.train_dir / "train.h5ad")
        model = celltypist.train(
            train.X,
            labels=train.obs["cell_type"].to_numpy(),
            genes=train.var_names.to_numpy(),
            use_SGD=USE_SGD,
            feature_selection=FEATURE_SELECTION,
            top_genes=TOP_GENES,
            balance_cell_type=BALANCE_CELL_TYPE,
            epochs=EPOCHS,
            n_jobs=N_JOBS,
            random_state=SEED,
        )
    model.write(str(args.out_dir / "model.pkl"))

    # Annotate: one call gives the per-cell and the majority-vote labels ---------------------
    result = celltypist.annotate(
        query, model=model, majority_voting=True, over_clustering=CLUSTER_COLUMN
    )
    frame = result.predicted_labels
    confidence = result.to_adata(
        insert_labels=False, insert_conf=True, insert_conf_by="predicted_labels"
    ).obs["conf_score"]
    cells = frame.assign(uncertainty=1.0 - confidence.reindex(frame.index))
    for name, column in {"per_cell": "predicted_labels", "majority": "majority_voting"}.items():
        folder = args.out_dir / name
        folder.mkdir()
        table = cells[[column, "uncertainty"]].rename(columns={column: "predicted_label"})
        table.rename_axis("barcode").reset_index().to_csv(
            folder / "cell_predictions.tsv", sep="\t", index=False
        )


if __name__ == "__main__":
    main()
