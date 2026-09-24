#!/bin/bash
# Etapa 4 / Tarea 2: BADF (Python, 1 core). Uso: sbatch run_badf.sbatch
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -J badf
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
/scratch/negishi/aaabt/venv_structural/bin/python badf.py
