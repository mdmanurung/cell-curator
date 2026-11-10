#!/bin/bash
# Timing run for protocol v7 cap: one whole annotation level, both marker methods.
#
# Run:     sbatch benchmarks/scripts/time_level_comparisons.sbatch RUN_DIR
# Inputs:  the whole-transcriptome AIDA slice (see time_level_comparisons.py)
# Outputs: RUN_DIR/steps.tsv
#
#SBATCH --job-name=cc-time-level
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
cd /exports/para-lipg-hpc/mdmanurung/cell-curator
.venv/bin/python benchmarks/scripts/time_level_comparisons.py "$1"
