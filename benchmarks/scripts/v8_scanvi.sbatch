#!/bin/bash
# scANVI comparator for protocol v8: train on the AIFI cells, map a query with scArches.
#
# Run:     sbatch benchmarks/scripts/v8_scanvi.sbatch TRAIN_DIR QUERY_H5AD OUT_DIR [TRUTH_TSV]
#          TRUTH_TSV given: also write the neutral table (a scoring step)
# Inputs:  TRAIN_DIR [v8_prepare_training.py], QUERY_H5AD input.h5ad
# Outputs: OUT_DIR/query.h5ad          query with a recovered counts layer
#          OUT_DIR/train/model/        trained scANVI reference
#          OUT_DIR/map/cell_predictions.tsv  barcode, predicted_label, uncertainty
#          with TRUTH_TSV OUT_DIR/neutral.tsv
# REUSE_TRAIN=OTHER_OUT_DIR/train (env var): link that finished reference instead of
# retraining; the training is seeded and does not see the query, so the model is the same.
# scANVI cannot abstain: the label is the arg-max soft probability; uncertainty = 1 - max.
#
#SBATCH --job-name=cc-v8-scanvi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1g.24gb:1
#SBATCH --exclude=res-hpc-gpu[01-02]
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=10:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
REPO=/exports/para-lipg-hpc/mdmanurung/cell-curator
ENV=/exports/para-lipg-hpc/mdmanurung/envs/cc-gpu
# v6 settings (scanvi.sbatch): 200 epochs for scVI and for scANVI, 100 for scArches, seed 0
TRAIN_EPOCHS=200
MAP_EPOCHS=100
SEED=0

nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
"$ENV/bin/python" -c "import scvi, torch; print('scvi', scvi.__version__, 'torch', torch.__version__, torch.version.cuda)"
mkdir "$3"
"$ENV/bin/python" "$REPO/benchmarks/scripts/v8_prepare_query.py" "$2" "$3/query.h5ad"
if [ -n "${REUSE_TRAIN:-}" ]; then
  ln -s "$(realpath "$REUSE_TRAIN")" "$3/train"
else
  "$ENV/bin/cell-curator" reference train-scanvi --reference "$1/train.h5ad" --labels-key cell_type \
    --counts-layer counts --max-epochs "$TRAIN_EPOCHS" --seed "$SEED" --output-dir "$3/train"
fi
"$ENV/bin/cell-curator" reference map-scarches --query "$3/query.h5ad" \
  --model-dir "$3/train/model" --max-epochs "$MAP_EPOCHS" --output-dir "$3/map"
if [ -n "${4:-}" ]; then
  "$REPO/.venv/bin/python" "$REPO/benchmarks/scripts/v8_to_neutral.py" \
    "$3/map/cell_predictions.tsv" "$2" "$4" "$1/label_map.tsv" "$3/neutral.tsv"
fi
