#!/bin/bash
# CellTypist comparator for protocol v8: retrain on the AIFI cells, annotate a query.
#
# Run:     sbatch benchmarks/scripts/v8_celltypist.sbatch TRAIN_DIR QUERY_H5AD OUT_DIR [TRUTH_TSV]
#          TRUTH_TSV given: also write the neutral tables (a scoring step; leave it out until due)
# Inputs:  see v8_celltypist.py
# Outputs: OUT_DIR/{model.pkl,per_cell,majority}/; with TRUTH_TSV, OUT_DIR/*/neutral.tsv
#
#SBATCH --job-name=cc-v8-celltypist
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
REPO=/exports/para-lipg-hpc/mdmanurung/cell-curator
ENV=/exports/para-lipg-hpc/mdmanurung/envs/cc-celltypist
"$ENV/bin/python" "$REPO/benchmarks/scripts/v8_celltypist.py" "$1" "$2" "$3"
if [ -n "${4:-}" ]; then
  for lane in per_cell majority; do
    "$REPO/.venv/bin/python" "$REPO/benchmarks/scripts/v8_to_neutral.py" \
      "$3/$lane/cell_predictions.tsv" "$2" "$4" "$1/label_map.tsv" "$3/$lane/neutral.tsv"
  done
fi
