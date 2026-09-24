#!/bin/bash
#SBATCH --job-name=qgrad_test
#SBATCH --output=slurm_test.%j.out
#SBATCH --error=slurm_test.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:20:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

for seg in 01 08 09 12; do
  echo "=== TEST seg${seg} ==="
  srun --mpi=pmix -n "$SLURM_NTASKS" lmp -in in.test_seg${seg}.lmp -log log.test_seg${seg}.${SLURM_JOB_ID}
  echo "=== fin TEST seg${seg} (exit $?) ==="
done
