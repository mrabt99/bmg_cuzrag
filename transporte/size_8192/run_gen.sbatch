#!/bin/bash
#SBATCH --job-name=s8192_gen
#SBATCH --output=gen.%j.out
#SBATCH --error=gen.%j.err
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
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.gen8192 -log log.gen.$SLURM_JOB_ID
