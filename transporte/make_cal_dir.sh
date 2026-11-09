#!/bin/bash
# Uso: ./make_cal_dir.sh <T> <tasa> <restart> <nsteps_npt> [nprod]
# Crea cal_<T>K_<tasa>/ con in.npt, in.nve, sbatch y scripts de analisis (copiados de cal_300K_2e10). NO lanza nada.
set -euo pipefail
T=$1; R=$2; RST=$3; NS=$4; NP=${5:-1000000}
D=cal_${T}K_$R; mkdir -p $D
cp cal_300K_2e10/{in.npt,in.nve,analisis_gk.py,analisis_replicas.py} $D/
sed -e "s|Vidrio 300 K, tasa 2e10: NPT 100 ps desde el vidrio relajado (10 ns a 300 K)|${T} K, tasa $R: NPT desde $RST|" \
    -e "s|gk300_npt_2e10|gk${T}_npt_$R|" \
    -e "s|-var rst [^ ]*|-var rst $RST|" -e "s|-var T 300.0|-var T ${T}.0|" -e "s|-var nsteps 100000|-var nsteps $NS|" \
    cal_300K_2e10/run_npt.sbatch > $D/run_npt.sbatch
sed -e "s|Vidrio 300 K, tasa 2e10|${T} K, tasa $R|" -e "s|gk300_nve_2e10|gk${T}_nve_$R|" \
    -e "s|-var T 300.0|-var T ${T}.0|" -e "s|-var nprod 1000000|-var nprod $NP|" \
    cal_300K_2e10/run_nve.sbatch > $D/run_nve.sbatch
echo "creado $D"; grep -h "srun" $D/run_npt.sbatch; grep -h "var T\|nprod" $D/run_nve.sbatch | tail -1
