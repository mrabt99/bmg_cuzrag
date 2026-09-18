#!/bin/bash
# =====================================================================
# Repeticion de la validacion v4 (Pdamp=10) con otra semilla de
# velocidades, para chequear si la explosion del primer intento
# (job 42846192, ~0.71 ns) es un evento raro estocastico o algo
# sistematico. Ver in.quench_1e10_validate_seed2.lmp y
# ESTADO_PROYECTO_meltquench_1e10.md.
# =====================================================================
#SBATCH --job-name=validate_1e10_seed2
#SBATCH --output=slurm_validate_seed2.%j.out
#SBATCH --error=slurm_validate_seed2.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=20:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.quench_1e10_validate_seed2.lmp -log log.validate_seed2.${SLURM_JOB_ID}
