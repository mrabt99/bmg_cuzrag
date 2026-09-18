#!/bin/bash
# =====================================================================
# Tramo 1/17 del templado a 1e10 K/s (v4: NPT, Pdamp=10.0), con la
# semilla de velocidades "buena" (seed2) y 2 ns extra de equilibracion
# a 2000 K antes de la rampa (sugerencia Jaroszewicz) -> 12 ns totales
# en este job en vez de 10. A la performance real medida en las
# validaciones (~27 ns/dia, 48 ranks) son ~10.7 h; incluso al ritmo
# conservador del log del profesor (6.845 ns/dia) serian ~42 h -> se
# pide 48 h de margen igual. Pensado para encadenarse con
# --dependency=afterok via chain_submit.sh (NO afterany: si este job
# muere a mitad de camino no queremos que el siguiente arranque
# leyendo un restart inexistente).
# =====================================================================
#SBATCH --job-name=quench1e10_seg01
#SBATCH --output=slurm_seg01.%j.out
#SBATCH --error=slurm_seg01.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=48:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.quench_1e10_seg01.lmp -log log.seg01.${SLURM_JOB_ID}
