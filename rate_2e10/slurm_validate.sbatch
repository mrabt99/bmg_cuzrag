#!/bin/bash
# =====================================================================
# Validacion del fix v4 (NPT, Pdamp=10.0) -- 4.5 ns, ANTES de lanzar
# la cadena de produccion. Ver in.quench_1e10_validate.lmp e
# ESTADO_PROYECTO_meltquench_1e10.md (v1 exploto a 1.67 ns con
# Pdamp=1.0; esto corre 2.7x mas alla de ese punto para confirmar que
# el fix aguanta).
#
# Estimado: 4.5e6 pasos / 79.2 pasos/s (perf medida en 48 ranks para
# este sistema) =~ 15.8 h -> se pide 20 h de margen.
# =====================================================================
#SBATCH --job-name=validate_1e10
#SBATCH --output=slurm_validate.%j.out
#SBATCH --error=slurm_validate.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=20:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.quench_1e10_validate.lmp -log log.validate.${SLURM_JOB_ID}
