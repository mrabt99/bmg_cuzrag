#!/bin/bash
#SBATCH -J ext_rylt
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONDONTWRITEBYTECODE=1
PY=/scratch/negishi/aaabt/venv_structural/bin/python
$PY -c "import scipy" 2>/dev/null || { module purge; module load anaconda/2024.02-py311; PY=python; }
$PY ryltsev_cmp.py
