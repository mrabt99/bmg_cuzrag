#!/bin/bash
# Descriptores de Laves para las configs nuevas de qgrad_ext (OVITO via anaconda) y despues el analisis. Escribe SOLO en qgrad_ext/laves/.
#SBATCH -J ext_laves
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
module purge; module load anaconda/2024.02-py311
export PYTHONDONTWRITEBYTECODE=1
E=/scratch/negishi/aaabt/bmg_cuzrag/qgrad_ext
for l in 14 15 16 17u 17r; do python laves_features.py qgrad_ext $l $E/snaps/$l.dump $E/boo/$l.boo.dump; done
python laves_analysis_ext.py
