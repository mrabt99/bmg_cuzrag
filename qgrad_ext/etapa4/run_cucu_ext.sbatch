#!/bin/bash
# Cu-Cu multi-frame (5 curvas) + profundidad del valle por T de la rama. Escribe SOLO en qgrad_ext/etapa4/.
#SBATCH -J ext_cucu
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONDONTWRITEBYTECODE=1
PY=/scratch/negishi/aaabt/venv_structural/bin/python
$PY -c "import scipy" 2>/dev/null || { module purge; module load anaconda/2024.02-py311; PY=python; }
$PY cucu_multiframe_ext.py
