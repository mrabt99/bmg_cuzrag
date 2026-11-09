#!/bin/bash
# Etapa 0 / Etapa 2 (BOO) para qgrad. Copia adaptada de batch_restart_to_dump.sbatch y batch_boo.sbatch
# (los originales no se tocan; --mpi=pmi2 porque pmix esta roto; QOS standby).
#
# Uso: sbatch --job-name=qg_dump batch_qgrad.sbatch dump <manifest>
#      sbatch --job-name=qg_boo  batch_qgrad.sbatch boo  <manifest>
#   <manifest>: lineas "etiqueta ruta_al_restart".
#   dump -> snaps/<etiqueta>.dump      (restart_to_dump.lmp: id type x y z, wrapped)
#   boo  -> boo/<etiqueta>.boo.dump    (validate_boo.lmp: Q4 Q6 W4hat W6hat para nnn 12/14/15/16)
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
MODE="$1"; MANIFEST="$2"
SA=/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis
POT=/scratch/negishi/aaabt/bmg_cuzrag/rate_2e10/ZrCuAg.lammps.eam
mkdir -p snaps boo logs
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
while read -r label rpath; do
    [ -z "$label" ] && continue
    if [ "$MODE" = dump ]; then in=$SA/restart_to_dump.lmp; out=snaps/${label}.dump
    else in=$SA/boo/validate_boo.lmp; out=boo/${label}.boo.dump; fi
    echo "### [$MODE] $label <- $rpath"
    srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in "$in" \
        -var restart_file "$rpath" -var pot_file "$POT" -var out_file "$out" \
        -log "logs/${MODE}.${label}.${SLURM_JOB_ID}.log" < /dev/null
done < "$MANIFEST"
echo "### Done ($MODE, $MANIFEST)"
