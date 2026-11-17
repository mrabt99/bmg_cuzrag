#!/bin/bash
# Fase 3 (otras tasas, 1e11 y 1e12): crea las carpetas cal_<T>K_<tasa>/. NO lanza nada (para eso: lanzar_fase3.sh).
# Plantilla: cal_1100K_2e10 (in.npt/in.nve con halts por dmin < 1.5 A y densidad < 6.0; analisis por trapecio).
# Restarts verificados 2026-11-10 con la ultima linea de thermo de cada segmento:
#   1e11 prod_1e11 segN termina a 1900-100N K: seg8 1105, seg9 1000, seg13 599, seg14 498, seg15 402 K;
#        glass300K = seg16 + 10 ns a 300 K (Pdamp 10), igual que 2e10.
#   1e12 rate_1e12_reference segN termina a 2000-100N K: seg9 1099, seg10 1007, seg14 598, seg15 496, seg16 401 K.
# NPT: 100 ps en el vidrio (300-600 K) y 200 ps en el liquido (1000-1100 K), igual que con 2e10.
# Uso (desde transporte/): ./preparar_fase3.sh      (saltea carpetas que ya existen)
set -euo pipefail
cd "$(dirname "$0")"
TPL=cal_1100K_2e10
R11=../../rate_1e11/prod_1e11
R12=../../rate_1e12_reference
#     T     tasa  restart                                         pasos NPT   tanda
CASOS="300  1e11  $R11/restart.glass300K_N65536_1e11              100000      A
       400  1e11  $R11/restart.quench_N65536_1e11.seg15           100000      A
       500  1e11  $R11/restart.quench_N65536_1e11.seg14           100000      A
       600  1e11  $R11/restart.quench_N65536_1e11.seg13           100000      A
       400  1e12  $R12/restart.quench_N65536_1e12.seg16           100000      A
       500  1e12  $R12/restart.quench_N65536_1e12.seg15           100000      A
       600  1e12  $R12/restart.quench_N65536_1e12.seg14           100000      A
       1000 1e11  $R11/restart.quench_N65536_1e11.seg9            200000      B
       1100 1e11  $R11/restart.quench_N65536_1e11.seg8            200000      B
       1000 1e12  $R12/restart.quench_N65536_1e12.seg10           200000      B
       1100 1e12  $R12/restart.quench_N65536_1e12.seg9            200000      B
       700  1e11  $R11/restart.quench_N65536_1e11.seg12           200000      C
       800  1e11  $R11/restart.quench_N65536_1e11.seg11           200000      C
       900  1e11  $R11/restart.quench_N65536_1e11.seg10           200000      C
       700  1e12  $R12/restart.quench_N65536_1e12.seg13           200000      C
       800  1e12  $R12/restart.quench_N65536_1e12.seg12           200000      C
       900  1e12  $R12/restart.quench_N65536_1e12.seg11           200000      C"
echo "$CASOS" | while read -r T R RST NS TANDA; do
  D=cal_${T}K_$R
  if [ -e "$D" ]; then echo "$D: ya existe -> salteado"; continue; fi
  # RST es relativo a la carpeta nueva (../../); desde transporte/ es un ../ menos
  if [ ! -f "${RST#../}" ]; then echo "$D: NO existe el restart $RST -> salteado"; continue; fi
  mkdir "$D"
  cp $TPL/{in.npt,in.nve,analisis_gk.py,analisis_replicas.py} "$D"/
  sed -e "s|^# 1100 K, tasa 2e10: NPT desde [^ ]*\.|# ${T} K, tasa $R (Fase 3, tanda $TANDA): NPT desde $RST.|" \
      -e "s|gk1100_npt_2e10|gk${T}_npt_$R|" \
      -e "s|-var rst [^ ]*|-var rst $RST|" -e "s|-var T 1100.0|-var T ${T}.0|" -e "s|-var nsteps 200000|-var nsteps $NS|" \
      $TPL/run_npt.sbatch > "$D"/run_npt.sbatch
  sed -e "s|^# 1100 K, tasa 2e10:|# ${T} K, tasa $R (Fase 3, tanda $TANDA):|" -e "s|gk1100_nve_2e10|gk${T}_nve_$R|" \
      -e "s|-var T 1100.0|-var T ${T}.0|" \
      $TPL/run_nve.sbatch > "$D"/run_nve.sbatch
  echo "creado $D (tanda $TANDA)"
  grep -h "^srun" "$D"/run_npt.sbatch; grep -h "job-name" "$D"/run_npt.sbatch "$D"/run_nve.sbatch; grep -h -- "-var T" "$D"/run_nve.sbatch
done
