#!/bin/bash
# Fase 2 (eta vs T a 2e10): NPT 200 ps + 3 replicas NVE 1 ns (semillas 1001/2002/3003) en 1100, 1300-1600 K.
# Restarts: prod_2e10 segN termina a (1900-100N) K -> 1100=seg8, 1300=seg6, 1400=seg5, 1500=seg4, 1600=seg3.
# in.npt/in.nve llevan halts por dmin < 1.5 A y densidad < 6.0 (ademas del de T).
# Uso (desde transporte/): ./lanzar_fase2_2e10.sh [T ...]   (sin argumentos: las 5 T). Lanzar SOLO con go-ahead.
set -euo pipefail
cd "$(dirname "$0")"
TS=${*:-"1100 1300 1400 1500 1600"}
for T in $TS; do
  d=cal_${T}K_2e10
  # proteccion contra relanzamientos (2026-11-06 se lanzo 3 veces): saltear si ya hay salidas o jobs en cola
  if ls $d/cal_r*_hflux.dat >/dev/null 2>&1 || squeue -u $USER -h -o %j | grep -q "^gk${T}_n[pv]t_2e10$"; then
    echo "$d: ya tiene salidas o jobs en cola -> salteado (borrar/mover las salidas a mano si se quiere rehacer)"; continue
  fi
  ( cd $d
    j=$(sbatch --parsable run_npt.sbatch)
    ids=$j
    for r in "1 1001" "2 2002" "3 3003"; do set -- $r
      ids="$ids $(sbatch --parsable --dependency=afterok:$j --export=ALL,REP=$1,SEED=$2 run_nve.sbatch)"; done
    echo "$d: NPT + NVE r1-r3 -> $ids" )
done
