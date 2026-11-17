#!/bin/bash
# VDOS del vidrio 2e10 a 300 K: NVT 20 ps + NVE 80 ps (20 orígenes de VACF de 4 ps). ~6-10 min.
# Usa la misma configuración NPT que cal_300K_2e10 (restart + V medio).
# Uso: sbatch run_vdos.sbatch
#SBATCH --job-name=vdos300_2e10
#SBATCH --output=vdos.%j.out
#SBATCH --error=vdos.%j.err
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
CAL=../cal_300K_2e10
VMEAN=$(awk '{print $4}' ${CAL}/cal_vmean.inc)
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.vdos -log log.vdos.$SLURM_JOB_ID \
  -var tag v300_2e10 -var seed 4004 -var T 300.0 -var vmean ${VMEAN} -var rst ${CAL}/cal_npt.restart
