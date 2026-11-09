#!/bin/bash
# Tasa graduada hacia gamma_min -- RAMA ext, tramo 09 (700->600 K, 1.500e+09 K/s).
# 66666667 pasos =~ 61.7h a ~300 pasos/s (48 ranks) -> se pide
# 5-16:19:00 de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=qgradext_seg09
#SBATCH --output=slurm_seg09.%j.out
#SBATCH --error=slurm_seg09.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=5-16:19:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_ext_seg09.lmp -log log.seg09.${SLURM_JOB_ID}
