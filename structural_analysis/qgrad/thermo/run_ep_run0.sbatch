#!/bin/bash
# Ep/atomo (run 0) de los 4 vidrios relajados a 300 K: 2e10, 1e11, 1e12, qgrad. Standby, 8 ranks, segundos.
#SBATCH -J qg_ep0
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:20:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
B=/scratch/negishi/aaabt/bmg_cuzrag
module purge; module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
mkdir -p ep_run0_logs
for t in 2e10:$B/rate_2e10/prod_2e10/restart.glass300K_N65536_2e10 1e11:$B/rate_1e11/prod_1e11/restart.glass300K_N65536_1e11 \
         1e12:$B/rate_1e12_reference/restart.glass300K_N65536_1e12 qgrad:$B/structural_analysis/qgrad/relax300K/restart.glass300K_N65536_qgrad; do
    tag=${t%%:*}; rp=${t#*:}
    srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in ep_run0.lmp -var restart_file "$rp" -var tag "$tag" -log ep_run0_logs/$tag.log < /dev/null
done
grep -h "^EP_RUN0" ep_run0_logs/*.log | sort > ep_run0_resultados.txt
cat ep_run0_resultados.txt
