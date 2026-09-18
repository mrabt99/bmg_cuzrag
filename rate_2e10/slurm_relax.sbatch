#!/bin/bash
# =====================================================================
# Etapa 3: relajacion final 300 K / 0 Pa, 10 ns (1e7 pasos, ~35 h ->
# se pide 48 h de margen). Se lanza DESPUES de seg17 con
# --dependency=afterok en chain_submit.sh.
# =====================================================================
#SBATCH --job-name=relax300K_1e10
#SBATCH --output=slurm_relax.%j.out
#SBATCH --error=slurm_relax.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=48:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.relax300K_1e10 -log log.relax.${SLURM_JOB_ID}
