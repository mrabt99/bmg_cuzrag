#!/bin/bash
# Verificacion fix ave/correlate vs numpy: NVT 10 ps + NVE 0.1 ns, 65536 atomos (~10 min).
#SBATCH --job-name=verif_acf
#SBATCH --output=verif_acf.%j.out
#SBATCH --error=verif_acf.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
source_vmean=$(awk '{print $4}' ../cal_1200K_2e10/cal_vmean.inc)
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.verif -log log.verif.$SLURM_JOB_ID -var vmean ${source_vmean}
