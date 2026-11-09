#!/bin/bash
# Ensayo de traccion PILOTO: 1e12, 300 K, 1e9 1/s, x, hasta strain 0.2 (250000 pasos).
# Estimacion: ~313 pasos/s con 48 ranks (medido en la relajacion 2e10) -> ~13 min;
# con deform + dumps ~15-20 min. Walltime 1 h. QOS standby (jobs < 4 h).
# Header copiado de qgrad_ext/slurm_relax.sbatch (corrio bien), con --mpi=pmi2.
#SBATCH --job-name=tens_pilot
#SBATCH --output=slurm_tens_pilot.%j.out
#SBATCH --error=slurm_tens_pilot.%j.err
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
SAMPLE=1e12; TEMP=300; ERATE=1e9
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.tension_pilot \
     -var sample $SAMPLE -var T $TEMP -var erate $ERATE \
     -log log.tens_${SAMPLE}_${TEMP}K_${ERATE}_x.${SLURM_JOB_ID}
