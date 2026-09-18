#!/bin/bash
# Produccion 2e10 K/s -- tramo 9/16. 5e6 pasos =~ 4.6h a ~300 pasos/s
# (48 ranks) -> se pide 12h de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=q2e10_seg09
#SBATCH --output=slurm_seg09.%j.out
#SBATCH --error=slurm_seg09.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=12:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.quench_2e10_seg09.lmp -log log.seg09.${SLURM_JOB_ID}
