#!/bin/bash
# popV comparator for protocol v5: train on donor split A, annotate split B.
#
# Run:     sbatch benchmarks/scripts/popv.sbatch WORK_DIR ENV_DIR REPO_DIR [HVG]
# Inputs:  WORK_DIR/splits/{A,B}.h5ad   [export_splits.py]
#          ENV_DIR   venv with popv, scvi-tools, and numpy<2
#          REPO_DIR  checkout at the protocol commit, for run_popv.py
# Outputs: WORK_DIR/popv/cell_predictions.tsv
#
#SBATCH --job-name=cc-v5-popv
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1g.24gb:1
#SBATCH --exclude=res-hpc-gpu[01-02]
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/cc-v5-work/logs/%x-%j.out

set -euo pipefail
nvidia-smi --query-gpu=name,driver_version --format=csv
"$2/bin/python" -c "from importlib.metadata import version as v; print('popv', v('popv'), 'scvi-tools', v('scvi-tools'))"
"$2/bin/python" "$3/benchmarks/scripts/run_popv.py" "$1/splits" "$1/popv" ${4:-}
