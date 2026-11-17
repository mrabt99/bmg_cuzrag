#!/bin/bash
# VDOS de muy baja frecuencia: el mismo vidrio 2e10 a 300 K replicado 2x2x2 (524 288 átomos, L = 209,5 Å).
# Primer modo transversal de la caja: 0,077 THz (contra 0,154 THz sin replicar).
# NVT 20 ps + NVE 400 ps = 10 orígenes de VACF de 40 ps (t_c hasta ~38 ps, resolución ~0,02 THz).
# Costo estimado: 23 min (65 536 át., 48 núcleos, 420 ps) x 8 átomos x 48/128 núcleos ≈ 70 min.
# Salidas: v300_2e10_rep2_vacf_<i>.dat, v300_2e10_rep2_thermo.dat. Análisis:
#   analiza_reducida.py v300_2e10_rep2 5 10 20 30 38
# Uso: sbatch run_vdos_rep2.sbatch
#SBATCH --job-name=vdosR2_300_2e10
#SBATCH --output=vdosR2.%j.out
#SBATCH --error=vdosR2.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=03:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
CAL=../cal_300K_2e10
VMEAN=$(awk '{print $4}' ${CAL}/cal_vmean.inc)
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.vdos -log log.vdosR2.$SLURM_JOB_ID \
  -var tag v300_2e10_rep2 -var seed 6006 -var T 300.0 -var vmean ${VMEAN} -var rst ${CAL}/cal_npt.restart \
  -var rep 2 -var LSEG 40000 -var NORIG 10 -var sv 4
