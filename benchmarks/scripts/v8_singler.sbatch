#!/bin/bash
# SingleR comparator for protocol v8: train on the AIFI cells, classify a query.
#
# Run:     sbatch benchmarks/scripts/v8_singler.sbatch TRAIN_DIR QUERY_H5AD OUT_DIR [TRUTH_TSV]
#          TRUTH_TSV given: also write the neutral table (a scoring step)
# Inputs:  TRAIN_DIR [v8_prepare_training.py], QUERY_H5AD input.h5ad
# Outputs: OUT_DIR/{train_export,query_export}/, OUT_DIR/cell_predictions.tsv,
#          with TRUTH_TSV OUT_DIR/neutral.tsv
#
#SBATCH --job-name=cc-v8-singler
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
REPO=/exports/para-lipg-hpc/mdmanurung/cell-curator
PY=/exports/para-lipg-hpc/mdmanurung/envs/cc-celltypist/bin/python
export PATH=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin:$PATH
# SingleR 2.12.0 sits in the user library; envs/R-cc-singler is empty
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5

mkdir "$3"
"$PY" "$REPO/benchmarks/scripts/v8_export_for_r.py" "$1/train.h5ad" "$3/train_export"
"$PY" "$REPO/benchmarks/scripts/v8_export_for_r.py" "$2" "$3/query_export"
Rscript "$REPO/benchmarks/scripts/v8_singler.R" "$3/train_export" "$3/query_export" "$3/cell_predictions.tsv"
if [ -n "${4:-}" ]; then
  "$REPO/.venv/bin/python" "$REPO/benchmarks/scripts/v8_to_neutral.py" \
    "$3/cell_predictions.tsv" "$2" "$4" "$1/label_map.tsv" "$3/neutral.tsv"
fi
