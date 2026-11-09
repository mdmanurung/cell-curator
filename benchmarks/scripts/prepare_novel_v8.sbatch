#!/bin/bash
# Protocol v8: cluster one novel dataset (HVG, scale, PCA, kNN, Leiden) for the annotators.
#
# Run:     sbatch benchmarks/scripts/prepare_novel_v8.sbatch DATASET RUN_DIR
# Inputs:  see prepare_novel_v8.py
# Outputs: RUN_DIR/ (see prepare_novel_v8.py)
#
#SBATCH --job-name=cc-prep-v8
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
cd /exports/para-lipg-hpc/mdmanurung/cell-curator
/exports/para-lipg-hpc/mdmanurung/envs/cc-v8/bin/python benchmarks/scripts/prepare_novel_v8.py "$1" "$2"
