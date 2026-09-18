#!/bin/bash
# =====================================================================
# Tramo 11/17 del templado a 1e10 K/s (v4: NPT, Pdamp=10.0).
# 1e7 pasos (10 ns) =~ 35 h a 79.2 pasos/s medidos (48 ranks) -> se
# pide 48 h de margen. Pensado para encadenarse con
# --dependency=afterok via chain_submit.sh (NO afterany: si este job
# muere a mitad de camino no queremos que el siguiente arranque
# leyendo un restart inexistente).
# =====================================================================
#SBATCH --job-name=quench1e10_seg11
#SBATCH --output=slurm_seg11.%j.out
#SBATCH --error=slurm_seg11.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=48:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.quench_1e10_seg11.lmp -log log.seg11.${SLURM_JOB_ID}
