#!/bin/bash
# popV comparator for protocol v8: train on the AIFI cells, annotate a query.
#
# Run:     sbatch benchmarks/scripts/v8_popv.sbatch TRAIN_DIR QUERY_H5AD OUT_DIR [TRUTH_TSV]
#          TRUTH_TSV given: also write the neutral table (a scoring step)
# Inputs:  TRAIN_DIR [v8_prepare_training.py], QUERY_H5AD input.h5ad
# Outputs: OUT_DIR/splits/{A,B}.h5ad   reference (A, a link to train.h5ad) and query (B)
#          OUT_DIR/popv/cell_predictions.tsv  barcode, predicted_label, uncertainty
#          with TRUTH_TSV OUT_DIR/neutral.tsv
# Reuses run_popv.py unchanged (v6 method set, batch key donor_id, uncertainty =
# 1 - agreeing methods / methods run). popV has no abstention rule: every cell it keeps
# gets its most-voted label. A cell popV drops in preprocessing is Unknown in the neutral
# table (--barcodes), as in v6 (coverage 0.971).
#
#SBATCH --job-name=cc-v8-popv
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1g.24gb:1
#SBATCH --exclude=res-hpc-gpu[01-02]
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
REPO=/exports/para-lipg-hpc/mdmanurung/cell-curator
ENV=/exports/para-lipg-hpc/mdmanurung/envs/cc-popv
HVG=4000 # v6: popV's default 4,000 highly variable genes

nvidia-smi --query-gpu=name,driver_version --format=csv
"$ENV/bin/python" -c "from importlib.metadata import version as v; print('popv', v('popv'), 'scvi-tools', v('scvi-tools'))"
mkdir "$3" "$3/splits"
ln -s "$(realpath "$1/train.h5ad")" "$3/splits/A.h5ad"
"$ENV/bin/python" "$REPO/benchmarks/scripts/v8_prepare_query.py" "$2" "$3/splits/B.h5ad"
"$ENV/bin/python" "$REPO/benchmarks/scripts/run_popv.py" "$3/splits" "$3/popv" "$HVG"
if [ -n "${4:-}" ]; then
  "$REPO/.venv/bin/python" "$REPO/benchmarks/scripts/v8_to_neutral.py" \
    "$3/popv/cell_predictions.tsv" "$2" "$4" "$1/label_map.tsv" "$3/neutral.tsv"
fi
