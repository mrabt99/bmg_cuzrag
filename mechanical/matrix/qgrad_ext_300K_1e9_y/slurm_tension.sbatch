#!/bin/bash
# Traccion: muestra qgrad_ext, T = 300 K, erate = 1e9 s^-1, eje y. Pasos estimados: 250000 (~14 min a 300 pasos/s). standby (<4 h), pmi2, 48 ranks.
#SBATCH --job-name=t_qgrad_ext_300_1e9_y
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.tension -var sample qgrad_ext -var T 300 -var erate 1e9 -var dir y \
     -log log.tens_qgrad_ext_300K_1e9_y.${SLURM_JOB_ID}
