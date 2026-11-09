#!/bin/bash
# Etapa 4 (RDF) para qgrad: los 4 tipos (total, partial, sametype, cross) sobre los restarts de un manifest.
# Reusa SIN EDITAR los run_rdf*.lmp de etapa4/ (parametrizados por -var). --mpi=pmi2, standby.
# Uso: sbatch --job-name=qg_rdf batch_rdf_qgrad.sbatch <manifest>     (lineas "etiqueta ruta_restart")
# Salida: <tipo>/raw/qgrad/<etiqueta>.rdf.txt y .log
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=01:30:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
MANIFEST="$1"
E4=/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis/etapa4
POT=/scratch/negishi/aaabt/bmg_cuzrag/rate_2e10/ZrCuAg.lammps.eam
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
declare -A IN=( [rdf]=$E4/rdf/run_rdf.lmp [rdf_partial]=$E4/rdf_partial/run_rdf_partial.lmp
                [rdf_sametype]=$E4/rdf_sametype/run_rdf_sametype.lmp [rdf_cross]=$E4/rdf_cross/run_rdf_cross.lmp )
for kind in rdf rdf_partial rdf_sametype rdf_cross; do
    mkdir -p "$kind/raw/qgrad"
    while read -r label rpath; do
        [ -z "$label" ] && continue
        echo "### [$kind] $label <- $rpath"
        srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in "${IN[$kind]}" \
            -var restart_file "$rpath" -var pot_file "$POT" -var out_file "$kind/raw/qgrad/${label}.rdf.txt" \
            -log "$kind/raw/qgrad/${label}.log" < /dev/null
    done < "$MANIFEST"
done
echo "### Done RDF x4 ($MANIFEST)"
