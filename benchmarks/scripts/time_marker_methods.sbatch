#!/bin/bash
# Timing run for tracker row 19: CPU Wilcoxon and logistic regression on whole clusters.
#
# Run:     sbatch benchmarks/scripts/time_marker_methods.sbatch RUN_DIR
# Inputs:  the whole-transcriptome AIDA slice (see time_marker_methods.py)
# Outputs: RUN_DIR/timing.tsv, RUN_DIR/versions.txt
#
#SBATCH --job-name=cc-time-markers
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
cd /exports/para-lipg-hpc/mdmanurung/cell-curator
.venv/bin/python benchmarks/scripts/time_marker_methods.py "$1"
