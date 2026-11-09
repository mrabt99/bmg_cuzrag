#!/bin/bash
# Traccion: muestra 2e10, T = 900 K, erate = 1e9 s^-1, eje y. Pasos estimados: 850000 (~47 min a 300 pasos/s). standby (<4 h), pmi2, 48 ranks.
#SBATCH --job-name=t_2e10_900_1e9_y
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.tension -var sample 2e10 -var T 900 -var erate 1e9 -var dir y \
     -log log.tens_2e10_900K_1e9_y.${SLURM_JOB_ID}
