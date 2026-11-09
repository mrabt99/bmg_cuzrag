#!/bin/bash
# Voronoi default de OVITO sobre qgrad/snaps (necesita anaconda: venv_structural no tiene ovito).
#SBATCH -J qg_voro
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR/aggregate"
module purge; module load anaconda/2024.02-py311
export PYTHONDONTWRITEBYTECODE=1
python voronoi_qgrad.py "$@"
