#!/bin/bash
# PRUEBA CORTA elasticas 300 K: NPT 2000 pasos, NVT 500 + 1000 pasos por estado, solo columnas 1 (xx) y 4 (yz). Sin significado estadistico: solo valida el flujo.
#SBATCH --job-name=el300_smoke
#SBATCH --output=slurm_smoke.%j.out
#SBATCH --error=slurm_smoke.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:20:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
EXTRA="-var nnpt 2000 -var nequil 500 -var nrun 1000 -var nser 250" DIRS="1 4" bash run_all.sh
