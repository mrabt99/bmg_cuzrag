#!/bin/bash
# Tasa graduada hacia gamma_min -- tramo 7/12 (900->800K, 2.172e+09 K/s).
# 46047198 pasos =~ 42.6h a ~300 pasos/s (48 ranks) -> se pide
# 3-22:18:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgrad_seg07
#SBATCH --output=slurm_seg07.%j.out
#SBATCH --error=slurm_seg07.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=3-22:18:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_grad_seg07.lmp -log log.seg07.${SLURM_JOB_ID}
