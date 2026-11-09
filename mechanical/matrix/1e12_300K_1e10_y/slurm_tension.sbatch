#!/bin/bash
# Traccion: muestra 1e12, T = 300 K, erate = 1e10 s^-1, eje y. Pasos estimados: 70000 (~4 min a 300 pasos/s). standby (<4 h), pmi2, 48 ranks.
#SBATCH --job-name=t_1e12_300_1e10_y
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
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
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.tension -var sample 1e12 -var T 300 -var erate 1e10 -var dir y \
     -log log.tens_1e12_300K_1e10_y.${SLURM_JOB_ID}
