#!/bin/bash
# Prueba stress/atom vs centroid/stress/atom en heat/flux (vidrio 300 K, 2e10). NVT 2 ps + NVE 5 ps, ~1-2 min.
#SBATCH --job-name=verif_centroid
#SBATCH --output=verif_centroid.%j.out
#SBATCH --error=verif_centroid.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:20:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
vmean=$(awk '{print $4}' ../cal_300K_2e10/cal_vmean.inc)
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.centroid -log log.centroid.$SLURM_JOB_ID -var vmean ${vmean}
