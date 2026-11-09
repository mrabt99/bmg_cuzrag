#!/bin/bash
# PRUEBA CORTA elasticas 0 K: preparacion con tolerancias REALES (para medir iteraciones/tiempo) + columnas 1 (xx) y 4 (yz), up = 1e-4.
#SBATCH --job-name=el0K_smoke
#SBATCH --output=slurm_smoke.%j.out
#SBATCH --error=slurm_smoke.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:30:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
UPS="1e-4" DIRS="1 4" bash run_0K.sh
