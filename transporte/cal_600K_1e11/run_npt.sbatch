#!/bin/bash
# 600 K, tasa 1e11 (Fase 3, tanda A): NPT desde ../../rate_1e11/prod_1e11/restart.quench_N65536_1e11.seg13. Genera cal_npt.restart y cal_vmean.inc
#SBATCH --job-name=gk600_npt_1e11
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
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.npt -log log.npt.$SLURM_JOB_ID -var rst ../../rate_1e11/prod_1e11/restart.quench_N65536_1e11.seg13 -var tag cal -var T 600.0 -var nsteps 100000
