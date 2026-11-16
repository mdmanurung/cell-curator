#!/usr/bin/env bash
#SBATCH --partition=short
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=cc_sif
# Requested: short partition (1 h limit), 40 min, 4 CPUs, 8 GB (the last build used 6 min).
# Run from this folder: sbatch -o _build/sbatch_%j.log build_sif.sbatch [NAME]
module load container/apptainer/1.2.5/gcc-8.5.0
cd "$SLURM_SUBMIT_DIR"
bash build_sif.sh "${1:-consensus}"
