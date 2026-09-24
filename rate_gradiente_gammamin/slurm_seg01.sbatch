#!/bin/bash
# Tasa graduada hacia gamma_min -- tramo 1/12 (1500->1400K, 2.000e+10 K/s).
# 5000000 pasos =~ 4.6h a ~300 pasos/s (48 ranks) -> se pide
# 10:42:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgrad_seg01
#SBATCH --output=slurm_seg01.%j.out
#SBATCH --error=slurm_seg01.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=10:42:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.quench_grad_seg01.lmp -log log.seg01.${SLURM_JOB_ID}
