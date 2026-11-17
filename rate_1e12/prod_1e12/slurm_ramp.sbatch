#!/bin/bash
# Produccion 1e12 K/s UNIFICADA: los 16 tramos (1900 -> 300 K) uno tras otro en UN job.
# Cada tramo = 1e5 pasos (~6 min a ~290 pasos/s con 48 ranks) -> ~1.7 h en total; se piden 4 h (standby).
# Si un tramo falla (fix halt -> error), set -e corta el job; para retomar desde el tramo N:
#   sbatch --export=ALL,FROM=N slurm_ramp.sbatch
#SBATCH --job-name=q1e12_ramp
#SBATCH --output=slurm_ramp.%j.out
#SBATCH --error=slurm_ramp.%j.err
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

FROM=${FROM:-1}
for n in $(seq "$FROM" 16); do
  i=$(printf "%02d" "$n")
  echo "### $(date '+%F %T') tramo $i"
  srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.quench_1e12_seg${i}.lmp -log log.seg${i}.${SLURM_JOB_ID}
  test -s restart.quench_N65536_1e12.seg${n}
done
echo "### $(date '+%F %T') rampa completa"
