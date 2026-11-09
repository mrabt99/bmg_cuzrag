#!/bin/bash
# Elasticas a 300 K del vidrio 2e10. Pasos: NPT 50 ps + (NVT 20 + 100 ps) x (1 referencia + 12 perturbaciones) = ~1.6e6 pasos;
# a ~245 pasos/s (caja triclinica, medido en la prueba corta) son ~110 min -> walltime 3 h. standby (<4 h). pmi2.
#SBATCH --job-name=el300_1e12
#SBATCH --output=slurm_elastic300K.%j.out
#SBATCH --error=slurm_elastic300K.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=03:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
bash run_all.sh
