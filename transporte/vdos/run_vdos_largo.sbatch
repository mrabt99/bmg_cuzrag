#!/bin/bash
# VDOS de baja frecuencia (búsqueda del pico bosónico), vidrio 2e10 a 300 K.
# Mismo in.vdos que run_vdos.sbatch, pero 20 orígenes de 20 ps (resolución < 0,1 THz) en vez de 4 ps.
# NVT 20 ps + NVE 400 ps ≈ 22 min (el job 45840090 tardó 5 min en 100 ps). Semilla nueva (5005).
# Salidas: v300_2e10_L20_vacf_<i>.dat, v300_2e10_L20_thermo.dat (no pisa las de 4 ps).
# Uso: sbatch run_vdos_largo.sbatch
#SBATCH --job-name=vdosL20_300_2e10
#SBATCH --output=vdosL20.%j.out
#SBATCH --error=vdosL20.%j.err
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
CAL=../cal_300K_2e10
VMEAN=$(awk '{print $4}' ${CAL}/cal_vmean.inc)
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.vdos -log log.vdosL20.$SLURM_JOB_ID \
  -var tag v300_2e10_L20 -var seed 5005 -var T 300.0 -var vmean ${VMEAN} -var rst ${CAL}/cal_npt.restart \
  -var LSEG 20000 -var NORIG 20 -var sv 4
