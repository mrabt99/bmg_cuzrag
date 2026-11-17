#!/bin/bash
# Elasticas a 0 K. La PREPARACION (minimizar atomos y caja, ~30 min) se hace UNA sola vez en prep/ y todas las magnitudes de deformacion la reutilizan:
# carpeta up_<valor>/ por cada una, con enlaces a los .mod y a prep/restart.equil.
# Uso: UPS="1e-4 5e-4" bash run_0K.sh      (EXTRA="-var maxiter 200 ..." y DIRS="1 4" para pruebas)
set -euo pipefail
HERE=$(pwd); NR=${SLURM_NTASKS:-48}; LMP="srun --mpi=pmi2 -n $NR lmp"
mkdir -p prep; cd prep; ln -sf ../*.mod .
$LMP -in ../in.elastic_prep ${EXTRA:-} -log log.prep
cd "$HERE"
for u in ${UPS:-1e-4 5e-4}; do
    mkdir -p up_$u; cd up_$u; ln -sf ../*.mod .; ln -sf ../prep/restart.equil .
    for d in ${DIRS:-1 2 3 4 5 6}; do $LMP -in ../in.elastic_dir -var up $u -var dir $d ${EXTRA:-} -log log.dir$d; done
    cd "$HERE"
done
