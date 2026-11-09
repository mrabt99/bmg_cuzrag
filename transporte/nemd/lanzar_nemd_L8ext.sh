#!/bin/bash
# Extension de L = 419 A (n=8) con Q = 0.5 eV/ps: tramos 4-6 desde L8.part3.restart (+3 ns, total 5 ns de produccion).
# Objetivo: medir el ruido real de kappa a Q = 0.5 (superbloques de 500 ps) y decidir si los errores de
# analiza_nemd.py estan subestimados (ver plan, bloque 2026-11-09). ~1 h por tramo, encadenados con afterok.
# Uso (desde nemd/): ./lanzar_nemd_L8ext.sh   Lanzar SOLO con go-ahead.
set -euo pipefail
cd "$(dirname "$0")"
TAG=L8; P0=4; P1=6
[ -f ${TAG}.part3.restart ] || { echo "falta ${TAG}.part3.restart"; exit 1; }
if ls ${TAG}_p${P0}_Tprof.dat >/dev/null 2>&1 || squeue -u $USER -h -o %j | grep -q "^nemd300$"; then
  echo "$TAG: ya hay tramo $P0 o jobs nemd300 en cola -> no se lanza"; exit 1
fi
ids=""; dep=""
for p in $(seq $P0 $P1); do
  j=$(sbatch --parsable $dep --export=ALL,N=8,PART=$p,SEED=801,Q=0.5,TAG=$TAG run_nemd.sbatch)
  ids="$ids $j"; dep="--dependency=afterok:$j"
done
echo "$TAG: tramos $P0-$P1 ->$ids"
