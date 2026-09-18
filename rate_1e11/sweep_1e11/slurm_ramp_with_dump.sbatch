#!/bin/bash
# Regeneracion con dump del segmento 1 global (2000->1900K, 1e11 K/s).
# Mismos 16 tasks que in.ramp.lmp original para minimizar diferencias de
# punto flotante -- 1e6 pasos tardo 2:21:11 la primera vez, se pide 6h
# de margen. Job independiente, sin dependency.
#SBATCH --job-name=sweep_1e11_redump
#SBATCH --output=slurm_ramp_with_dump.%j.out
#SBATCH --error=slurm_ramp_with_dump.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=06:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.ramp_with_dump.lmp -log log.ramp_with_dump.${SLURM_JOB_ID}
