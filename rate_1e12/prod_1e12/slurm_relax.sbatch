#!/bin/bash
# Relajacion final 300 K / 10 ns (1e7 pasos), igual que 2e10 y 1e11 (la de 1e11 tardo 9 h 00 min con 48 ranks).
# Mas de 4 h -> QOS normal (sin standby). Va encadenado con --dependency=afterok al job de la rampa.
#SBATCH --job-name=relax300K_1e12
#SBATCH --output=slurm_relax.%j.out
#SBATCH --error=slurm_relax.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=14:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.relax300K_1e12 -log log.relax.${SLURM_JOB_ID}
