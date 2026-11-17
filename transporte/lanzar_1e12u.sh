#!/bin/bash
# Encola el transporte de la 1e12 unificada (carpetas de preparar_1e12u.sh): por T, NPT + 3 replicas NVE de 1 ns
# (semillas 1001/2002/3003, afterok del NPT), igual que lanzar_fase3.sh. Con tope: 3 carriles; dentro de un carril cada T
# arranca cuando termino la T anterior (afterany), asi corren a lo sumo 9 jobs de transporte a la vez. 40 jobs, ~40 nodo-h.
# Uso (desde transporte/): ./lanzar_1e12u.sh      Lanzar SOLO con go-ahead. Saltea carpetas con salidas o logs previos.
set -euo pipefail
cd "$(dirname "$0")"
CARRILES=("300 1200 600 900" "400 1100 700" "500 1000 800")
for c in "${CARRILES[@]}"; do
  prev=""
  for T in $c; do
    d=cal_${T}K_1e12u
    [ -d "$d" ] || { echo "$d: no existe -> salteado"; continue; }
    if compgen -G "$d/cal_npt.restart" >/dev/null || compgen -G "$d/cal_r*_hflux.dat" >/dev/null || compgen -G "$d/gk_n*.out" >/dev/null; then
      echo "$d: ya tiene salidas o logs -> salteado"; continue; fi
    dep=""; [ -n "$prev" ] && dep="--dependency=afterany:$prev"
    j=$(cd "$d" && sbatch --parsable $dep run_npt.sbatch); ids="$j"; nve=""
    for r in "1 1001" "2 2002" "3 3003"; do set -- $r
      k=$(cd "$d" && sbatch --parsable --dependency=afterok:$j --export=ALL,REP=$1,SEED=$2 run_nve.sbatch); ids="$ids $k"; nve="$nve:$k"; done
    prev="${nve#:}"
    echo "$d: NPT + NVE r1-r3 -> $ids"
  done
done
