#!/bin/bash
# Prueba de linealidad: NEMD a L = 419 A (n=8) con Q = 0.25 eV/ps (la mitad), tag L8q025, semilla 811.
# 5 tramos de 1 ns encadenados con afterok (1 ns de transitorio + 4 ns de produccion; ~1 h por tramo).
# Uso (desde nemd/): ./lanzar_nemd_q025.sh   Lanzar SOLO con go-ahead.
set -euo pipefail
cd "$(dirname "$0")"
TAG=L8q025; NPART=5
if ls ${TAG}_p1_Tprof.dat >/dev/null 2>&1 || squeue -u $USER -h -o %j | grep -q "^nemd300$"; then
  echo "$TAG: ya hay salidas o jobs nemd300 en cola -> no se lanza"; exit 1
fi
ids=""; dep=""
for p in $(seq 1 $NPART); do
  j=$(sbatch --parsable $dep --export=ALL,N=8,PART=$p,SEED=811,Q=0.25,TAG=$TAG run_nemd.sbatch)
  ids="$ids $j"; dep="--dependency=afterok:$j"
done
echo "$TAG: tramos 1-$NPART ->$ids"
