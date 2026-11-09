#!/bin/bash
# Test de rendimiento (NO es produccion): solo in.nve (2000 pasos NVT + 20000 NVE con heat/flux), 48 ranks.
# No incluye NPT: su velocidad ya se conoce de los tramos del enfriamiento.
#SBATCH --job-name=gk_perf
#SBATCH --output=gk_perf.%j.out
#SBATCH --error=gk_perf.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:30:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
# usa el restart original seg7 directamente (sin reescalar caja), prefijo "perf"
srun --mpi=pmi2 -n $SLURM_NTASKS lmp -in in.nve -log log.perf_nve.$SLURM_JOB_ID \
  -var tag perf -var nnvt 2000 -var nprod 20000 -var rst ../../rate_1e11/prod_1e11/restart.quench_N65536_1e11.seg7
