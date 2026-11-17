#!/bin/bash
# Fase 3: encola NPT + 3 replicas NVE 1 ns (semillas 1001/2002/3003, afterok) en las carpetas de preparar_fase3.sh.
#   tanda A: kappa en el vidrio  -> cal_{300,400,500,600}K_1e11, cal_{400,500,600}K_1e12   (28 jobs, ~23 nodo-h)
#   tanda B: liquido sobreenfriado -> cal_{1000,1100}K_{1e11,1e12}                          (16 jobs, ~13 nodo-h)
# Uso (desde transporte/): ./lanzar_fase3.sh A | B | <carpeta> ...      Lanzar SOLO con go-ahead.
# Proteccion contra relanzamientos: saltea la carpeta si ya tiene salidas, logs de jobs previos o jobs en cola.
set -euo pipefail
cd "$(dirname "$0")"
case "${1:-}" in
  A) DIRS="cal_300K_1e11 cal_400K_1e11 cal_500K_1e11 cal_600K_1e11 cal_400K_1e12 cal_500K_1e12 cal_600K_1e12" ;;
  B) DIRS="cal_1000K_1e11 cal_1100K_1e11 cal_1000K_1e12 cal_1100K_1e12" ;;
  C) DIRS="cal_700K_1e11 cal_800K_1e11 cal_900K_1e11 cal_700K_1e12 cal_800K_1e12 cal_900K_1e12" ;;   # kappa cerca de Tg (24 jobs)
  "") echo "uso: $0 A | B | C | <carpeta> ..."; exit 1 ;;
  *) DIRS="$*" ;;
esac
COLA=$(squeue -u "$USER" -h -o %j)
for d in $DIRS; do
  [ -d "$d" ] || { echo "$d: no existe (correr preparar_fase3.sh) -> salteado"; continue; }
  jn=$(grep -h "job-name" "$d"/run_npt.sbatch | sed 's/.*=//; s/_npt_/_n(pt|ve)_/')
  if compgen -G "$d/cal_npt.restart" >/dev/null || compgen -G "$d/cal_r*_hflux.dat" >/dev/null || compgen -G "$d/gk_n*.out" >/dev/null \
     || echo "$COLA" | grep -qE "^${jn}$"; then
    echo "$d: ya tiene salidas, logs o jobs en cola -> salteado (mover a mano si se quiere rehacer)"; continue
  fi
  ( cd "$d"
    j=$(sbatch --parsable run_npt.sbatch); ids=$j
    for r in "1 1001" "2 2002" "3 3003"; do set -- $r
      ids="$ids $(sbatch --parsable --dependency=afterok:$j --export=ALL,REP=$1,SEED=$2 run_nve.sbatch)"; done
    echo "$d: NPT + NVE r1-r3 -> $ids" )
done
