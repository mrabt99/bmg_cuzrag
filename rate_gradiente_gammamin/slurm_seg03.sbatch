#!/bin/bash
# Tasa graduada hacia gamma_min -- tramo 3/12 (1300->1200K, 9.542e+09 K/s).
# 10480478 pasos =~ 9.7h a ~300 pasos/s (48 ranks) -> se pide
# 21:51:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgrad_seg03
#SBATCH --output=slurm_seg03.%j.out
#SBATCH --error=slurm_seg03.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=21:51:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_grad_seg03.lmp -log log.seg03.${SLURM_JOB_ID}
