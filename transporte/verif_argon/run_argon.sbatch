#!/bin/bash
# Verificacion de la normalizacion de kappa con el ejemplo de la doc de LAMMPS (Ar solido, 70 K; la doc dice ~0.29 W/mK).
# 2 versiones (real = doc verbatim; metal = ruta de unidades de nuestro pipeline) x 2 ensembles (nvt = doc, nve) x 3 semillas. 256 atomos.
#SBATCH --job-name=verif_argon
#SBATCH --output=verif_argon.%j.out
#SBATCH --error=verif_argon.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=02:00:00
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
for ens in nvt nve; do for seed in 102486 20241 777; do
  for v in real metal; do
    srun --mpi=pmi2 -n1 --exact lmp -in in.ar_$v -log log.ar_${v}_${ens}_$seed -var seed $seed -var ens $ens -var tag ar_${v}_${ens}_$seed > out.ar_${v}_${ens}_$seed 2>&1 &
  done
  wait
done; done
