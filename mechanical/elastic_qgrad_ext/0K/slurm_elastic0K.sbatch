#!/bin/bash
# Elasticas a 0 K del vidrio 2e10: delta = 1e-4 y verificacion a 5e-4 (carpetas up_1e-4/ y up_5e-4/).
# Preparacion UNA vez (~30 min medidos en la prueba corta) + 12 minimizaciones cortas (segundos) por cada delta -> walltime 1.5 h (margen: la preparacion ya rozo 30 min). standby, pmi2.
#SBATCH --job-name=el0K_qgrad_ext
#SBATCH --output=slurm_elastic0K.%j.out
#SBATCH --error=slurm_elastic0K.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=01:30:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
UPS="1e-4 5e-4" bash run_0K.sh
