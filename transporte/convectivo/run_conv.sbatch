#!/bin/bash
# Prueba del termino convectivo: NVT 20 ps + NVE 1 ns (~1 h) de una replica, 2e10, a la T dada.
# Parte de ../cal_<T>K_2e10/cal_npt.restart y de su volumen medio (solo lectura). Escribe en conv_<T>K_2e10/.
# Uso (desde transporte/convectivo/): sbatch --export=ALL,T=900,REP=1,SEED=7007 --job-name=conv900_r1 run_conv.sbatch
#SBATCH --output=conv.%j.out
#SBATCH --error=conv.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=03:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${T:?falta T}"; REP=${REP:-1}; SEED=${SEED:-7007}
SRC=$(cd ../cal_${T}K_2e10 && pwd); D=conv_${T}K_2e10; mkdir -p $D; cd $D
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
source_vmean=$(awk '{print $4}' $SRC/cal_vmean.inc)
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in ../in.nve_conv -log log.conv_r${REP}.$SLURM_JOB_ID \
  -var tag conv_r${REP} -var seed ${SEED} -var T ${T}.0 -var vmean ${source_vmean} -var rst $SRC/cal_npt.restart -var nprod 1000000
