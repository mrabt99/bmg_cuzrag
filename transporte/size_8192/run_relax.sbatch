#!/bin/bash
#SBATCH --job-name=s8192_relax
#SBATCH --output=relax.%j.out
#SBATCH --error=relax.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=04:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.relax8192 -log log.relax.$SLURM_JOB_ID
