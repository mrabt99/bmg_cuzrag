#!/bin/bash
# Regenera analisis_replicas (trapecio) en cada cal_* y refresca el unico png viejo (cal_1200K_1e11/cal_r1_gk.png). No lanza jobs de Slurm.
PY=/scratch/negishi/aaabt/venv_structural/bin/python
for d in cal_300K_2e10 cal_300K_1e12 cal_400K_2e10 cal_500K_2e10 cal_600K_2e10 cal_1000K_2e10 cal_1200K_2e10 cal_1200K_1e11 cal_1200K_1e12; do
  T=$(echo $d | sed 's/cal_\([0-9]*\)K_.*/\1/'); V=$(awk '{print $4}' $d/cal_vmean.inc)
  ( cd $d && echo "# analisis_replicas.py (trapecio, $(date +%F)) T=$T V=$V" > analisis_replicas_trapecio.txt && $PY analisis_replicas.py $T $V cal_r1 cal_r2 cal_r3 >> analisis_replicas_trapecio.txt 2>&1 ) && echo "ok $d"
done
( cd cal_1200K_1e11 && $PY analisis_gk.py cal_r1 1200 $(awk '{print $4}' cal_vmean.inc) > cal_r1_gk_trapecio.out 2>&1 && echo "ok png 1200K_1e11 r1" )
echo FIN
