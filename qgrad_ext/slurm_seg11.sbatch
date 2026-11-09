#!/bin/bash
# Tasa graduada hacia gamma_min -- tramo 11/12 (500->400K, 1.000e+11 K/s).
# 1000000 pasos =~ 0.9h a ~300 pasos/s (48 ranks) -> se pide
# 02:33:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgradext_seg11
#SBATCH --output=slurm_seg11.%j.out
#SBATCH --error=slurm_seg11.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=02:33:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_ext_seg11.lmp -log log.seg11.${SLURM_JOB_ID}
