#!/bin/bash
# Encola la prueba del termino convectivo: 2e10 a 600, 900 y 1200 K, 3 replicas de 1 ns (semillas 7007/8008/9009), 9 jobs sin NPT.
# Uso (desde transporte/convectivo/): ./lanzar_conv.sh [T ...]     Lanzar SOLO con go-ahead.
# Proteccion contra relanzamientos: saltea la replica si ya tiene salida, log o job en cola con su nombre.
set -euo pipefail
cd "$(dirname "$0")"
TS=${*:-"600 900 1200"}
COLA=$(squeue -u "$USER" -h -o %j)
for T in $TS; do
  [ -f ../cal_${T}K_2e10/cal_npt.restart ] || { echo "falta ../cal_${T}K_2e10/cal_npt.restart -> salteado"; continue; }
  for r in "1 7007" "2 8008" "3 9009"; do set -- $r
    jn=conv${T}_r$1; d=conv_${T}K_2e10
    if compgen -G "$d/conv_r$1_conv.dat" >/dev/null || compgen -G "$d/log.conv_r$1.*" >/dev/null || echo "$COLA" | grep -qx "$jn"; then
      echo "$jn: ya tiene salida, log o job en cola -> salteado"; continue
    fi
    echo "$jn -> $(sbatch --parsable --export=ALL,T=$T,REP=$1,SEED=$2 --job-name=$jn run_conv.sbatch)"
  done
done
