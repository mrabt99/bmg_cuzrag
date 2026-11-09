#!/bin/bash
# Calibracion 1200 K: NVT 20 ps + NVE 1 ns (~1.3 h), una replica. Requiere cal_npt.restart y cal_vmean.inc.
# Uso: sbatch run_nve.sbatch            (replica 1, semilla 1001)
#      sbatch --export=ALL,REP=2,SEED=2002 run_nve.sbatch
#SBATCH --job-name=gk_nve
#SBATCH --output=gk_nve.%j.out
#SBATCH --error=gk_nve.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=03:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
REP=${REP:-1}; SEED=${SEED:-1001}
source_vmean=$(awk '{print $4}' cal_vmean.inc)
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.nve -log log.nve_r${REP}.$SLURM_JOB_ID \
  -var tag cal_r${REP} -var seed ${SEED} -var T 1200.0 -var vmean ${source_vmean} -var rst cal_npt.restart -var nprod 1000000
