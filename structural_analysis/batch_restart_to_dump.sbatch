#!/bin/bash
# Etapa 0 (batch): convierte los 17 restarts de una tasa en
# snaps/01.dump ... snaps/17.dump, via restart_to_dump.lmp.
#
# Uso: sbatch --job-name=dump_<tasa> batch_restart_to_dump.sbatch <tasa>
#   <tasa>/manifest.txt debe tener 17 lineas (paths absolutos a
#   restarts, en orden 01..17).
#SBATCH --output=slurm_batch_%x.%j.out
#SBATCH --error=slurm_batch_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=02:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

RATE="$1"
POT_FILE="/scratch/negishi/aaabt/bmg_cuzrag/rate_2e10/ZrCuAg.lammps.eam"
MANIFEST="${RATE}/manifest.txt"
SNAPS_DIR="${RATE}/snaps"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

mapfile -t restarts < "$MANIFEST"

i=0
for restart_path in "${restarts[@]}"; do
    i=$((i+1))
    n=$(printf "%02d" "$i")
    out_file="${SNAPS_DIR}/${n}.dump"
    echo "### [$RATE] ${n}.dump <- ${restart_path}"
    srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in restart_to_dump.lmp \
        -var restart_file "$restart_path" \
        -var pot_file "$POT_FILE" \
        -var out_file "$out_file" \
        -log "${RATE}/log.${n}.${SLURM_JOB_ID}" \
        < /dev/null
done

echo "### Done: $i restarts converted for $RATE"
