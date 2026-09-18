#!/bin/bash
# Regenerar el segmento 1 (2000->1900K) con dump de trayectoria, que
# faltaba en el sweep original. Mismos 16 tasks que el sweep original
# para reproducir la misma trayectoria lo mas fielmente posible.
#SBATCH --job-name=sweep_2e10_redump
#SBATCH --output=slurm_ramp_with_dump.%j.out
#SBATCH --error=slurm_ramp_with_dump.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=16:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.ramp_with_dump.lmp -log log.ramp_with_dump.${SLURM_JOB_ID}
