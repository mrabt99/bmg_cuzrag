#!/bin/bash
# Tablas y figura de la rama qgrad_ext (corre tras dump, BOO y Voronoi). Escribe SOLO en qgrad_ext/aggregate/.
#SBATCH -J ext_post
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR/aggregate"
export PYTHONDONTWRITEBYTECODE=1
PY=/scratch/negishi/aaabt/venv_structural/bin/python
$PY consolidate_ext.py
$PY regimen_ext.py
$PY fig_ext.py
