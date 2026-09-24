#!/bin/bash
# Etapa 4 / Tarea 4: RDF mismo-tipo (total + Zr-Zr + Cu-Cu) de los 17 restarts de UNA tasa.
# Uso: sbatch --job-name=rdf_sametype_<tasa> batch_rdf_sametype.sbatch <tasa>
#   lee ../manifests/<tasa>.txt (17 paths, orden 01..17)
#   escribe raw/<tasa>/NN.rdf.txt y raw/<tasa>/NN.log
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:30:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
RATE="$1"
POT_FILE="/scratch/negishi/aaabt/bmg_cuzrag/rate_2e10/ZrCuAg.lammps.eam"
MANIFEST="../manifests/${RATE}.txt"
mkdir -p "raw/${RATE}"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

mapfile -t restarts < "$MANIFEST"
i=0
for restart_path in "${restarts[@]}"; do
    i=$((i+1)); n=$(printf "%02d" "$i")
    echo "### [$RATE] ${n} <- ${restart_path}"
    srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in run_rdf_sametype.lmp \
        -var restart_file "$restart_path" -var pot_file "$POT_FILE" \
        -var out_file "raw/${RATE}/${n}.rdf.txt" \
        -log "raw/${RATE}/${n}.log" < /dev/null
done
echo "### Done: $i configs for $RATE"
