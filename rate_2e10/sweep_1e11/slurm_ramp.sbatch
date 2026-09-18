#!/bin/bash
# Sweep de tasa -- 1e11 K/s, tramo 2000->1900K, 1000000 pasos.
# Job independiente, SIN dependency (todo el sweep corre en paralelo).
#SBATCH --job-name=sweep_1e11
#SBATCH --output=slurm_ramp.%j.out
#SBATCH --error=slurm_ramp.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=08:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.ramp.lmp -log log.ramp.${SLURM_JOB_ID}
