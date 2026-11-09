#!/bin/bash
# Regenera tablas y figuras 4 curvas con 300 K = 17r (relajado). Corre tras dump/BOO/RDF/Voronoi/features de 17r.
#SBATCH -J qg_post17r
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONDONTWRITEBYTECODE=1
SA=/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis
PY=/scratch/negishi/aaabt/venv_structural/bin/python
module purge; module load anaconda/2024.02-py311
(cd etapa4 && python setup_4curvas.py)
(cd etapa4/rdf && $PY build_rdf_csv_and_figs_4curvas.py)
for k in partial sametype cross; do (cd etapa4/rdf_$k && $PY build_rdf_${k}_csv_and_figs_4curvas.py); done
(cd etapa4/rdf_sametype && $PY cucu_split_4curvas.py)
(cd etapa4/fig4_style && $PY build_fig4_style_4curvas.py)
(cd etapa4/badf && $PY badf_4curvas.py)
(cd aggregate && $PY consolidate_qgrad.py && $PY regimen_qgrad.py && python step12_fig_compare_rates_4curvas.py && python step13_fig_compare_rates_log_4curvas.py)
(cd laves && $PY laves_analysis.py && $PY laves_clusters_z16ico.py && $PY laves_gr.py && $PY laves_sq.py && $PY export_flags.py qgrad 17r $SA/qgrad/snaps/17r.dump)
