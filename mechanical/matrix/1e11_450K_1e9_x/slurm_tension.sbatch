#!/bin/bash
# Traccion: muestra 1e11, T = 450 K, erate = 1e9 s^-1, eje x. Pasos estimados: 400000 (~22 min a 300 pasos/s). standby (<4 h), pmi2, 48 ranks.
#SBATCH --job-name=t_1e11_450_1e9_x
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
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
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.tension -var sample 1e11 -var T 450 -var erate 1e9 -var dir x \
     -log log.tens_1e11_450K_1e9_x.${SLURM_JOB_ID}
