#!/bin/bash
# NEMD (fix ehex) a 300 K, un tramo de 1 ns. Variables: N (replicas en z), PART, SEED. Ver in.nemd y README.md.
# Uso: sbatch --export=ALL,N=8,PART=1,SEED=801 run_nemd.sbatch
#      opcionales: Q (eV/ps, default 0.5) y TAG (default L<N>); p.ej. N=8,Q=0.25,TAG=L8q025
#SBATCH --job-name=nemd300
#SBATCH --output=nemd.%j.out
#SBATCH --error=nemd.%j.err
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
N=${N:-8}; PART=${PART:-1}; SEED=${SEED:-801}; Q=${Q:-0.5}; TAG=${TAG:-L${N}}
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.nemd -log log.${TAG}_p${PART}.$SLURM_JOB_ID \
  -var n ${N} -var tag ${TAG} -var part ${PART} -var seed ${SEED} -var Q ${Q}
