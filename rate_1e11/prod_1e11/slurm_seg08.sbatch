#!/bin/bash
# Produccion 1e11 K/s -- tramo 8/16. 1e6 pasos =~ 1h a ~300 pasos/s
# (48 ranks) -> se pide 4h de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=q1e11_seg08
#SBATCH --output=slurm_seg08.%j.out
#SBATCH --error=slurm_seg08.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=04:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.quench_1e11_seg08.lmp -log log.seg08.${SLURM_JOB_ID}
