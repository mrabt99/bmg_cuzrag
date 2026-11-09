#!/bin/bash
# BADF 4 curvas (Python, 1 core). Uso: sbatch run_badf_4curvas.sbatch
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -J qg_badf
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONDONTWRITEBYTECODE=1
/scratch/negishi/aaabt/venv_structural/bin/python badf_4curvas.py
