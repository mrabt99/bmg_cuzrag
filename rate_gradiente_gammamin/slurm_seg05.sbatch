#!/bin/bash
# Tasa graduada hacia gamma_min -- tramo 5/12 (1100->1000K, 4.552e+09 K/s).
# 21968082 pasos =~ 20.3h a ~300 pasos/s (48 ranks) -> se pide
# 1-21:15:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgrad_seg05
#SBATCH --output=slurm_seg05.%j.out
#SBATCH --error=slurm_seg05.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=1-21:15:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_grad_seg05.lmp -log log.seg05.${SLURM_JOB_ID}
