#!/bin/bash
# Regla 1: md5 de TODOS los archivos de structural_analysis/ salvo qgrad/.
# Uso: sbatch checksums.sbatch antes    (o: despues)
# Escribe checksums_<antes|despues>.txt (ordenado por ruta, para poder hacer diff).
#SBATCH --job-name=qgrad_md5
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
WHEN="${1:?uso: sbatch checksums.sbatch antes|despues}"
cd ..
find . -type f ! -path './qgrad/*' -print0 | sort -z | xargs -0 md5sum > "qgrad/checksums_${WHEN}.txt.tmp"
mv "qgrad/checksums_${WHEN}.txt.tmp" "qgrad/checksums_${WHEN}.txt"
echo "### $(wc -l < qgrad/checksums_${WHEN}.txt) archivos -> qgrad/checksums_${WHEN}.txt"
