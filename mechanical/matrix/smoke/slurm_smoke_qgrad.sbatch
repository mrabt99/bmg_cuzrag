#!/bin/bash
# PRUEBA CORTA (2e10 como regresion, qgrad y qgrad_ext como nuevas): x, 300 K, 1e9, nhold 1000, emax 0.005 (~1 min cada una)
#SBATCH --job-name=tens_qgrad_smoke
#SBATCH --output=slurm_qgrad_smoke.%j.out
#SBATCH --error=slurm_qgrad_smoke.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:20:00
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
for s in 2e10 qgrad qgrad_ext; do
    mkdir -p ${s}_x; cd ${s}_x
    srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in ../../in.tension -var sample $s -var T 300 -var erate 1e9 -var dir x -var nhold 1000 -var emax 0.005 -log log.${s}.${SLURM_JOB_ID} || echo "### FALLO EN $s"
    cd ..
done
