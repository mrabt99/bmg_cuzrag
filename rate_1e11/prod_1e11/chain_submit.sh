#!/bin/bash
# Encadena los 16 tramos de produccion (1900->300K a 1e11 K/s) + la
# relajacion final con --dependency=afterok (no afterany).
# Uso:
#   ./chain_submit.sh                # lanza seg01..seg16 + relax
#   ./chain_submit.sh --from 5       # arranca en seg05
set -euo pipefail
cd "$(dirname "$0")"

FROM=1
if [ "${1:-}" = "--from" ]; then
  FROM="$2"
fi

PREV=""
for i in $(seq -w "$FROM" 16); do
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
echo "el .err/.out del job que fallo y relanzar con: ./chain_submit.sh --from N"
