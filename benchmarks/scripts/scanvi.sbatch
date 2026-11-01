#!/bin/bash
# scANVI comparator for protocol v5: train on donor split A, map split B with scArches.
#
# Run:     sbatch benchmarks/scripts/scanvi.sbatch WORK_DIR ENV_DIR
# Inputs:  WORK_DIR/splits/{A,B}.h5ad              [export_splits.py]
#          ENV_DIR   venv with cell-curator at the protocol commit, scvi-tools, cupy
# Outputs: WORK_DIR/scanvi/train/model/            trained scANVI reference
#          WORK_DIR/scanvi/map/cell_predictions.tsv  barcode, predicted_label, uncertainty
#
#SBATCH --job-name=cc-v5-scanvi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1g.24gb:1
#SBATCH --exclude=res-hpc-gpu[01-02]
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
WORK=$1
ENV=$2
CLI="$ENV/bin/cell-curator"

nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
"$ENV/bin/python" -c "import scvi, torch; print('scvi', scvi.__version__, 'torch', torch.__version__, torch.version.cuda)"

"$CLI" reference train-scanvi --reference "$WORK/splits/A.h5ad" --labels-key cell_type \
  --counts-layer counts --max-epochs 200 --seed 0 --output-dir "$WORK/scanvi/train"
"$CLI" reference map-scarches --query "$WORK/splits/B.h5ad" \
  --model-dir "$WORK/scanvi/train/model" --max-epochs 100 --output-dir "$WORK/scanvi/map"
