#!/bin/bash
# Encola las cadenas NEMD: por cada n (replicas en z) encadena NPART tramos de 1 ns con afterok.
# n=4 (L=210 A, 32768 at.): 2 tramos; n=8 (419 A, 65536): 3 tramos; n=16 (838 A, 131072): 5 tramos.
# Uso (desde nemd/): ./lanzar_nemd.sh [n ...]   (sin argumentos: 4 8 16). Lanzar SOLO con go-ahead.
set -euo pipefail
cd "$(dirname "$0")"
declare -A NPART=([4]=2 [8]=3 [16]=5)
for n in ${*:-4 8 16}; do
  if ls L${n}_p1_Tprof.dat >/dev/null 2>&1 || squeue -u $USER -h -o "%j %k" | grep -q nemd300; then
    echo "n=$n: ya hay salidas o jobs nemd300 en cola -> salteado"; continue
  fi
  ids=""; dep=""
  for p in $(seq 1 ${NPART[$n]}); do
    j=$(sbatch --parsable $dep --export=ALL,N=$n,PART=$p,SEED=$((100*n+1)) run_nemd.sbatch)
    ids="$ids $j"; dep="--dependency=afterok:$j"
  done
  echo "n=$n: tramos 1-${NPART[$n]} ->$ids"
done
