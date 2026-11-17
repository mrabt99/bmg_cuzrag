#!/bin/bash
# Transporte de la 1e12 UNIFICADA (rate_1e12/, etiqueta 1e12u; 17/11/2026): crea cal_<T>K_1e12u/ para las MISMAS 10 T que
# tenia la 1e12 vieja (300-1200 K). NO lanza nada (para eso: lanzar_1e12u.sh). Copia de la logica de preparar_fase3.sh.
# Plantilla: cal_1100K_2e10 (halts por dmin < 1.5 A y densidad < 6.0; analisis por trapecio).
# Restarts: en rate_1e12/prod_1e12 segN termina a 1900-100N K (como 1e11; verificado 17/11 con la ultima linea de
# cada log: seg7 1199, seg8 1102, seg9 1002, seg10 900, seg11 801, seg12 700, seg13 601, seg14 501, seg15 401 K).
# NPT: 100 ps en el vidrio (300-600 K) y 200 ps desde 700 K, igual que la 1e12 vieja y 2e10/1e11.
set -euo pipefail
cd "$(dirname "$0")"
TPL=cal_1100K_2e10
R=../../rate_1e12/prod_1e12
#     T     restart                                       pasos NPT
CASOS="300  $R/restart.glass300K_N65536_1e12              100000
       400  $R/restart.quench_N65536_1e12.seg15           100000
       500  $R/restart.quench_N65536_1e12.seg14           100000
       600  $R/restart.quench_N65536_1e12.seg13           100000
       700  $R/restart.quench_N65536_1e12.seg12           200000
       800  $R/restart.quench_N65536_1e12.seg11           200000
       900  $R/restart.quench_N65536_1e12.seg10           200000
       1000 $R/restart.quench_N65536_1e12.seg9            200000
       1100 $R/restart.quench_N65536_1e12.seg8            200000
       1200 $R/restart.quench_N65536_1e12.seg7            200000"
echo "$CASOS" | while read -r T RST NS; do
  D=cal_${T}K_1e12u
  if [ -e "$D" ]; then echo "$D: ya existe -> salteado"; continue; fi
  if [ ! -f "${RST#../}" ]; then echo "$D: NO existe el restart $RST -> salteado"; continue; fi
  mkdir "$D"
  cp $TPL/{in.npt,in.nve,analisis_gk.py,analisis_replicas.py} "$D"/
  sed -e "s|^# 1100 K, tasa 2e10: NPT desde [^ ]*\.|# ${T} K, tasa 1e12u (1e12 unificada): NPT desde $RST.|" \
      -e "s|gk1100_npt_2e10|gk${T}_npt_1e12u|" \
      -e "s|-var rst [^ ]*|-var rst $RST|" -e "s|-var T 1100.0|-var T ${T}.0|" -e "s|-var nsteps 200000|-var nsteps $NS|" \
      $TPL/run_npt.sbatch > "$D"/run_npt.sbatch
  sed -e "s|^# 1100 K, tasa 2e10:|# ${T} K, tasa 1e12u (1e12 unificada):|" -e "s|gk1100_nve_2e10|gk${T}_nve_1e12u|" \
      -e "s|-var T 1100.0|-var T ${T}.0|" \
      $TPL/run_nve.sbatch > "$D"/run_nve.sbatch
  echo "creado $D: $(grep -o -- '-var rst [^ ]*' $D/run_npt.sbatch) $(grep -o -- '-var T [0-9.]*' $D/run_npt.sbatch) $(grep -o -- '-var nsteps [0-9]*' $D/run_npt.sbatch) | NVE: $(grep -o -- '-var T [0-9.]*' $D/run_nve.sbatch)"
done
