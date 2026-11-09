#!/bin/bash
#SBATCH -J qg_feat17r
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
module purge; module load anaconda/2024.02-py311
export PYTHONDONTWRITEBYTECODE=1
SA=/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis
python laves_features.py qgrad 17r $SA/qgrad/snaps/17r.dump $SA/qgrad/boo/17r.boo.dump
