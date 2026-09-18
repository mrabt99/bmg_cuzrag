#!/bin/bash
# Tarea E - validacion: corre el patron de referencia del director
# (~/legacy_scripts/boo_seba_ejemplo.in, solo Q4/Q6 nnn=12) y la
# version extendida (validate_boo.lmp, 4 nnn + wl/hat) sobre EL MISMO
# restart (2e10, seg8), para comparar columna a columna antes de
# escalar a los 51 restarts.
#SBATCH --output=slurm_validate.%j.out
#SBATCH --error=slurm_validate.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

echo "### baseline (patron del director, solo Q4/Q6 nnn=12)"
srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in /scratch/negishi/aaabt/legacy_scripts/boo_seba_ejemplo.in \
    -log baseline_seba.log.${SLURM_JOB_ID} < /dev/null

echo "### extendido (4 nnn + wl/hat)"
srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in validate_boo.lmp \
    -log validate.log.${SLURM_JOB_ID} < /dev/null

echo "### Done"
