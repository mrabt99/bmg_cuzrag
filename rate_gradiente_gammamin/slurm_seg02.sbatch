#!/bin/bash
# Tasa graduada hacia gamma_min -- tramo 2/12 (1400->1300K, 1.381e+10 K/s).
# 7238949 pasos =~ 6.7h a ~300 pasos/s (48 ranks) -> se pide
# 15:15:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgrad_seg02
#SBATCH --output=slurm_seg02.%j.out
#SBATCH --error=slurm_seg02.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=15:15:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_grad_seg02.lmp -log log.seg02.${SLURM_JOB_ID}
