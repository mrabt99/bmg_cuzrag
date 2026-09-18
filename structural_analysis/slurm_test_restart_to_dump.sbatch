#!/bin/bash
# Caso de prueba de Etapa 0 (restart_to_dump): UN restart de 2e10, sin
# batch. Solo lee el restart y escribe un dump (write_dump, sin run),
# deberia terminar en minutos.
#SBATCH --job-name=test_restart2dump
#SBATCH --output=slurm_test.%j.out
#SBATCH --error=slurm_test.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in restart_to_dump.lmp -log log.test.${SLURM_JOB_ID}
