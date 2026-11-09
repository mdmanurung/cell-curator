#!/bin/bash
# Protocol v8: learn the corpus proposal from the Human Immune Health Atlas.
#
# Run:     sbatch benchmarks/scripts/learn_corpus_v8.sbatch RUN_DIR   (RUN_DIR must not exist)
# Inputs:  see learn_corpus_v8.py
# Outputs: RUN_DIR/ (see learn_corpus_v8.py)
#
#SBATCH --job-name=cc-learn-v8
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=06:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
cd /exports/para-lipg-hpc/mdmanurung/cell-curator
.venv/bin/python benchmarks/scripts/learn_corpus_v8.py "$1"
