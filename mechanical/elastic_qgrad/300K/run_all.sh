#!/bin/bash
# Corre la etapa previa y despues las 6 columnas de Cij, una tras otra (cada una con 2 signos). Se lanza desde el sbatch.
# Variables opcionales:  EXTRA="-var nrun 1000 ..."  (se pasa a todos los lmp)   DIRS="1 4"  (solo esas columnas; para pruebas)
set -euo pipefail
NR=${SLURM_NTASKS:-48}
LMP="srun --mpi=pmi2 -n $NR lmp"
$LMP -in in.elastic_prep ${EXTRA:-} -log log.prep
for d in ${DIRS:-1 2 3 4 5 6}; do
    $LMP -in in.elastic_dir -var dir $d ${EXTRA:-} -log log.dir$d
done
