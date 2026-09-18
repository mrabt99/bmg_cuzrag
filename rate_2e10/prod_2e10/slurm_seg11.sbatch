#!/bin/bash
# Produccion 2e10 K/s -- tramo 11/16. 5e6 pasos =~ 4.6h a ~300 pasos/s
# (48 ranks) -> se pide 12h de margen. Encadenado con --dependency=afterok.
#SBATCH --job-name=q2e10_seg11
#SBATCH --output=slurm_seg11.%j.out
#SBATCH --error=slurm_seg11.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=12:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.quench_2e10_seg11.lmp -log log.seg11.${SLURM_JOB_ID}
