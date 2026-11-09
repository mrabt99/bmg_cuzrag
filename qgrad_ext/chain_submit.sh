#!/bin/bash
# Encadena los tramos de la rama qgrad_ext (plan A) con --dependency=afterok:
#   09: 700->600 K @1.5e9 K/s | 10a: 600->550 K @1.5e9 | 10: 550->500 K @1.5e9 (T_switch) | 11: 500->400 K @1e11 | 12: 400->300 K @1e11 | relax: 10 ns NPT a 300 K
# Copia adaptada de rate_gradiente_gammamin/chain_submit.sh. Si un tramo falla, la cadena se corta sola (afterok).
# Uso:
#   ./chain_submit.sh                       # lanza 09 10a 10 11 12 relax
#   ./chain_submit.sh --from 10a            # arranca en el tramo 10a (ej. tras un fallo)
#   ./chain_submit.sh --from 09 --upto 09   # solo un tramo
# IMPORTANTE: no se ejecuta solo; lo corre el usuario despues de revisar los inputs. La vigilancia de cristalizacion y de Ep(T) se hace
# con monitor_ext.py entre tramos (la decision de cortar la tasa lenta antes de 500 K es del usuario: hay restart en 550 K = fin de seg10a).
set -euo pipefail
cd "$(dirname "$0")"
STAGES=(09 10a 10 11 12 relax)
FROM=09; UPTO=relax
while [ $# -gt 0 ]; do
  case "$1" in
    --from) FROM="$2"; shift 2 ;;
    --upto) UPTO="$2"; shift 2 ;;
    *) echo "arg desconocido: $1"; exit 1 ;;
  esac
done
PREV=""; ON=0
for st in "${STAGES[@]}"; do
  [ "$st" = "$FROM" ] && ON=1
  [ "$ON" = 1 ] || continue
  if [ "$st" = relax ]; then script="slurm_relax.sbatch"; else script="slurm_seg${st}.sbatch"; fi
  [ -f "$script" ] || { echo "falta $script"; exit 1; }
  if [ -z "$PREV" ]; then jid=$(sbatch --parsable "$script"); else jid=$(sbatch --parsable --dependency=afterok:"$PREV" "$script"); fi
  echo "tramo $st -> job $jid  ($script)"
  PREV="$jid"
  [ "$st" = "$UPTO" ] && break
done
echo ""; echo "Cadena enviada. Seguir con: squeue -u \$USER"
echo "Recuperacion ante explosion: SIN hold isotermico (ver comentario en cada in.quench_ext_segNN.lmp)."
