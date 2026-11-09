#!/bin/bash
# NVE largo por tramos de 1 ns para eta a 800 K (ver in.nve_long). Variables: REP, SEED, PART.
#SBATCH --job-name=gk800_long
#SBATCH --output=gk_long.%j.out
#SBATCH --error=gk_long.%j.err
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
REP=${REP:-4}; SEED=${SEED:-4004}; PART=${PART:-1}
source_vmean=$(awk '{print $4}' cal_vmean.inc)
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.nve_long -log log.long_r${REP}_p${PART}.$SLURM_JOB_ID \
  -var tag long_r${REP} -var part ${PART} -var seed ${SEED} -var T 800.0 -var vmean ${source_vmean}
