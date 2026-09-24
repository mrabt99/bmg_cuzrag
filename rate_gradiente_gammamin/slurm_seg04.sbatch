#!/bin/bash
# Tasa graduada hacia gamma_min -- tramo 4/12 (1200->1100K, 6.590e+09 K/s).
# 15173529 pasos =~ 14.0h a ~300 pasos/s (48 ranks) -> se pide
# 1-07:25:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgrad_seg04
#SBATCH --output=slurm_seg04.%j.out
#SBATCH --error=slurm_seg04.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=1-07:25:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_grad_seg04.lmp -log log.seg04.${SLURM_JOB_ID}
