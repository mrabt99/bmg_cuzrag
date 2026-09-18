#!/bin/bash
# Tarea E (escalado): corre validate_boo.lmp (4 computes orientorder/atom,
# nnn 12/14/15/16, degrees 2 4 6, wl/hat yes) sobre los 17 restarts de UNA
# tasa, leidos de <rate>/manifest.txt (mismo archivo que ya usaba
# batch_restart_to_dump.sbatch para snaps/), potencial real EAM (no
# pair_style zero, no rerun -- validado contra el patron del director en
# boo_seba_ejemplo.in, ver validate_2e10_seg8.dump).
#
# Uso: sbatch --job-name=boo_<tasa> batch_boo.sbatch <tasa> <potencial>
#   <tasa>/manifest.txt debe tener 17 lineas (paths absolutos a restarts,
#   en orden 01..17, igual que para snaps/).
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
POT_FILE="$2"
MANIFEST="../${RATE}/manifest.txt"
OUT_DIR="${RATE}"

mkdir -p "$OUT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

mapfile -t restarts < "$MANIFEST"

i=0
for restart_path in "${restarts[@]}"; do
    i=$((i+1))
    n=$(printf "%02d" "$i")
    out_file="${OUT_DIR}/${n}.boo.dump"
    echo "### [$RATE] ${n}.boo.dump <- ${restart_path}"
    srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in validate_boo.lmp \
        -var restart_file "$restart_path" \
        -var pot_file "$POT_FILE" \
        -var out_file "$out_file" \
        -log "${OUT_DIR}/log.${n}.${SLURM_JOB_ID}" \
        < /dev/null
done

echo "### Done: $i restarts converted to BOO for $RATE"
