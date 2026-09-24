#!/bin/bash
# Encadena los 12 tramos (ventana graduada 1500->700K, tramos 1-8, +
# cola rapida 700->300K a 1e11 K/s, tramos 9-12) con --dependency=afterok.
# Uso:
#   ./chain_submit.sh                # lanza seg01..seg12
#   ./chain_submit.sh --from 5       # arranca en seg05
#   ./chain_submit.sh --from 5 --upto 8   # lanza solo seg05..seg08
#                                          # (util para pausar y chequear
#                                          # cristalizacion antes de seguir
#                                          # a los tramos mas lentos)
set -euo pipefail
cd "$(dirname "$0")"

FROM=1
UPTO=12
while [ $# -gt 0 ]; do
  case "$1" in
    --from) FROM="$2"; shift 2 ;;
    --upto) UPTO="$2"; shift 2 ;;
    *) echo "arg desconocido: $1"; exit 1 ;;
  esac
done

PREV=""
for i in $(seq -w "$FROM" "$UPTO"); do
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

echo ""
echo "Cadena enviada (tramos $FROM-$UPTO). Seguir con: squeue -u \$USER"
echo "Si un tramo falla, la cadena se corta sola (afterok); investigar"
echo "el .err/.out del job que fallo. Recuperacion en este proyecto:"
echo "SIN hold isotermico (ver comentario en cada in.quench_grad_segNN.lmp)."
echo ""
echo "IMPORTANTE (vigilancia de cristalizacion, tramos 1-8): revisar"
echo "check_crystallization.py sobre cada log.dat.segNN recien terminado"
echo "ANTES de lanzar el tramo mas lento siguiente -- no encadenar todo"
echo "1-8 a ciegas de una sola vez."
