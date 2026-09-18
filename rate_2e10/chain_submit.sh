#!/bin/bash
# =====================================================================
# Encadena los 17 tramos de templado + la relajacion final con
# --dependency=afterok (NO afterany: si un tramo muere a mitad de
# camino, el siguiente NO debe arrancar leyendo un restart inexistente
# o del tramo anterior -- lección de v3/sol67, ver
# ESTADO_PROYECTO_meltquench_1e10.md).
#
# Uso:
#   ./chain_submit.sh                # lanza seg01..seg17 + relax
#   ./chain_submit.sh --from 5        # arranca en seg05 (ya corrio 1-4)
#
# Antes de correr esto: confirmar que in.quench_1e10_validate.lmp
# termino SIN error (ni fix halt_T ni halt_P dispararon) via
# slurm_validate.sbatch. Si la validacion no corrio o fallo, no lanzar
# la cadena.
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

FROM=1
if [ "${1:-}" = "--from" ]; then
  FROM="$2"
fi

declare -A JOBID
PREV=""

for i in $(seq -w "$FROM" 17); do
  n=$((10#$i))
  script="slurm_seg${i}.sbatch"
  if [ -z "$PREV" ]; then
    jid=$(sbatch --parsable "$script")
  else
    jid=$(sbatch --parsable --dependency=afterok:"$PREV" "$script")
  fi
  echo "tramo $n -> job $jid  (input ${script})"
  PREV="$jid"
done

if [ -n "$PREV" ]; then
  jid=$(sbatch --parsable --dependency=afterok:"$PREV" slurm_relax.sbatch)
  echo "relax  -> job $jid  (input slurm_relax.sbatch)"
fi

echo ""
echo "Cadena enviada. Seguir con: squeue -u \$USER"
echo "Si un tramo falla, la cadena se corta sola (afterok); investigar"
echo "el .err del job que fallo y relanzar con: ./chain_submit.sh --from N"
