#!/bin/bash
# Descriptores de Laves para: referencias, qgrad (06..16,17u), controles 2e10 (01..17), 1e11 y 1e12 (todas las T).
#SBATCH -J qg_feat
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge; module load anaconda/2024.02-py311
export PYTHONDONTWRITEBYTECODE=1
SA=/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis
for t in C15 C14; do for k in 0 3 6 9; do python laves_features.py refs ${t}_300K_f$k refs/${t}_300K.pos.dump refs/${t}_300K.boo.dump $k; done
    python laves_features.py refs ${t}_0K refs/${t}_0K.pos.dump refs/${t}_0K.boo.dump 0; done
for l in 06 07 08 09 10 11 12 13 14 15 16 17u; do python laves_features.py qgrad $l $SA/qgrad/snaps/$l.dump $SA/qgrad/boo/$l.boo.dump; done
for r in 2e10 1e11 1e12; do for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17; do
    python laves_features.py $r $n $SA/$r/snaps/$n.dump $SA/boo/$r/$n.boo.dump; done; done
