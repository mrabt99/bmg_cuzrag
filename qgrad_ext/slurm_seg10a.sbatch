#!/bin/bash
# Tasa graduada hacia gamma_min -- RAMA ext, tramo 10a (600->550 K, 1.500e+09 K/s).
# 33333333 pasos =~ 30.9h a ~300 pasos/s (48 ranks) -> se pide
# 2-19:54:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgradext_seg10a
#SBATCH --output=slurm_seg10a.%j.out
#SBATCH --error=slurm_seg10a.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=2-19:54:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_ext_seg10a.lmp -log log.seg10a.${SLURM_JOB_ID}
