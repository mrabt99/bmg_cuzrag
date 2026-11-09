#!/bin/bash
# Vidrio 300 K, tasa 2e10: NPT 100 ps desde el vidrio relajado (10 ns a 300 K). Genera cal_npt.restart y cal_vmean.inc
#SBATCH --job-name=gk300_npt_2e10
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
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.npt -log log.npt.$SLURM_JOB_ID -var rst ../../rate_2e10/prod_2e10/restart.glass300K_N65536_2e10 -var tag cal -var T 300.0 -var nsteps 100000
