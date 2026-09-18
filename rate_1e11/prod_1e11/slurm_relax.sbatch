#!/bin/bash
# Relajacion final 300K/10ns -- cadena de produccion 1e11 K/s.
# 1e7 pasos, misma fisica/duracion que prod_2e10 (=~9.3h a ~300 pasos/s
# con 48 ranks) -> se pide 24h de margen.
#SBATCH --job-name=relax300K_1e11
#SBATCH --output=slurm_relax.%j.out
#SBATCH --error=slurm_relax.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=24:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.relax300K_1e11 -log log.relax.${SLURM_JOB_ID}
