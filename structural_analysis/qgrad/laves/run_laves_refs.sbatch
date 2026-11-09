#!/bin/bash
# Relaja y termaliza C15 y C14 de Cu2Zr con el EAM real (referencias del analisis de Laves).
#SBATCH -J qg_lavref
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge; module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
srun --mpi=pmi2 -n 4 lmp -in laves_ref.lmp -var tag C15 -var data refs/C15_init.data -var baro iso   -log refs/log.C15 < /dev/null
srun --mpi=pmi2 -n 4 lmp -in laves_ref.lmp -var tag C14 -var data refs/C14_init.data -var baro aniso -log refs/log.C14 < /dev/null
grep "^###" refs/log.C15 refs/log.C14
