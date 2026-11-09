#!/bin/bash
# 1000 K, tasa 1e12_N8192: NPT desde ../size_8192/restart.T1000_N8192_1e12. Genera cal_npt.restart y cal_vmean.inc
#SBATCH --job-name=gk1000_npt_1e12_N8192
#SBATCH --output=gk_npt.%j.out
#SBATCH --error=gk_npt.%j.err
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
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.npt -log log.npt.$SLURM_JOB_ID -var rst ../size_8192/restart.T1000_N8192_1e12 -var tag cal -var T 1000.0 -var nsteps 200000
