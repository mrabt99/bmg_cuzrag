#!/bin/bash
#SBATCH --job-name=qgrad_pmi2test
#SBATCH --output=slurm_test_pmi2.%j.out
#SBATCH --error=slurm_test_pmi2.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:15:00

cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623

echo "=== A: srun --mpi=pmi2 ==="
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.test_seg01.lmp -log log.A_pmi2.${SLURM_JOB_ID}
echo "=== fin A (exit $?) ==="

echo "=== B: mpirun ==="
mpirun -np "$SLURM_NTASKS" lmp -in in.test_seg01.lmp -log log.B_mpirun.${SLURM_JOB_ID}
echo "=== fin B (exit $?) ==="
