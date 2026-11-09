#!/bin/bash
# PRUEBAS CORTAS de in.tension (cada caso en su subcarpeta; el input es ../in.tension):
#   equiv_x_300K : x, 300 K, nhold 50000 (como el piloto), emax 0.005  -> comparar L0 y las primeras lineas de la curva con el piloto 2e10
#   dir_x/y/z_300K : nhold 1000, emax 0.005 en las tres direcciones
#   T450_x : x, 450 K con la rampa REAL (150000 pasos a 1e12 K/s), nhold 1000, emax 0.005
# Estimado: ~3 + 3x0.5 + ~10 min = ~15 min. standby, pmi2.
#SBATCH --job-name=tens_matrix_smoke
#SBATCH --output=slurm_matrix_smoke.%j.out
#SBATCH --error=slurm_matrix_smoke.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:40:00
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
run_case () {   # nombre T dir nhold emax
    mkdir -p "$1"; cd "$1"
    srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in ../../in.tension -var sample 2e10 -var T "$2" -var erate 1e9 -var dir "$3" \
         -var nhold "$4" -var emax "$5" -log "log.$1.${SLURM_JOB_ID}" || echo "### FALLO EN $1"
    cd ..
}
run_case equiv_x_300K  300 x 50000 0.005
run_case dir_x_300K    300 x 1000  0.005
run_case dir_y_300K    300 y 1000  0.005
run_case dir_z_300K    300 z 1000  0.005
run_case T450_x        450 x 1000  0.005
