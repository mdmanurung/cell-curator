#!/bin/bash
# Protocol v8: download the AIFI and Tabula Sapiens files (about 66 GB) and write the
# pre-run tables.
#
# Run:     sbatch benchmarks/scripts/acquire_v8.sbatch RUN_DIR   (RUN_DIR must not exist)
# Inputs:  CELLxGENE collections e9360edf-... and e5f58829-... (see acquire_v8.py)
# Outputs: RUN_DIR/ (see acquire_v8.py)
#
#SBATCH --job-name=cc-acquire-v8
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
cd /exports/para-lipg-hpc/mdmanurung/cell-curator
.venv/bin/python benchmarks/scripts/acquire_v8.py "$1"
