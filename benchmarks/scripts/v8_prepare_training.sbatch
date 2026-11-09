#!/bin/bash
# Training set for the v8 comparators: the sampled AIFI cells, all genes.
#
# Run:     sbatch benchmarks/scripts/v8_prepare_training.sbatch SAMPLED_TSV TRAIN_DIR
# Inputs:  see v8_prepare_training.py
# Outputs: TRAIN_DIR/train.h5ad, TRAIN_DIR/label_map.tsv
#
#SBATCH --job-name=cc-v8-train-prep
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=03:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
ENV=/exports/para-lipg-hpc/mdmanurung/envs/cc-celltypist
"$ENV/bin/python" /exports/para-lipg-hpc/mdmanurung/cell-curator/benchmarks/scripts/v8_prepare_training.py "$1" "$2"
