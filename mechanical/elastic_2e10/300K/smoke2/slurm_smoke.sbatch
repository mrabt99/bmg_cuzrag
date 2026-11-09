#!/bin/bash
# PRUEBA CORTA 2 (tras corregir ave running): conteos A PROPOSITO no multiplos de nser, para probar que el promedio ya no depende de la alineacion.
#SBATCH --job-name=el300_smoke2
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
EXTRA="-var nnpt 3000 -var nequil 700 -var nrun 3000 -var nser 1000" DIRS="1 4" bash run_all.sh
