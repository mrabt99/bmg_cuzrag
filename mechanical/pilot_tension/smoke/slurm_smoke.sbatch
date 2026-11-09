#!/bin/bash
# PRUEBA CORTA del input de traccion: nhold=1000, emax=0.005 -> 1000 + 5000 = 6000 pasos (~1 min).
# Mismo input (../in.tension_pilot); las salidas quedan en esta carpeta (smoke/).
#SBATCH --job-name=tens_smoke
#SBATCH --output=slurm_tens_smoke.%j.out
#SBATCH --error=slurm_tens_smoke.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:10:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
SAMPLE=2e10; TEMP=300; ERATE=1e9
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in ../in.tension_pilot \
     -var sample $SAMPLE -var T $TEMP -var erate $ERATE -var nhold 1000 -var emax 0.005 \
     -log log.tens_${SAMPLE}_${TEMP}K_${ERATE}_x.${SLURM_JOB_ID}
