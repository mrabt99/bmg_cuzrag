#!/bin/bash
# Relajacion 300K/10ns de qgrad. 1e7 pasos =~ 9.3h a ~300 pasos/s (48 ranks) -> 24h de margen.
# Mismo criterio que rate_2e10/prod_2e10/slurm_relax.sbatch, pero con --mpi=pmi2 (pmix roto).
#SBATCH --job-name=relax300K_qgrad
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
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.relax300K_qgrad -log log.relax.${SLURM_JOB_ID}
