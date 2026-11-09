#!/bin/bash
# Rehace SOLO el punto de 300 K de qgrad con el vidrio RELAJADO (17r) cuando termine el job de relajacion (relax300K/), y regenera
# tablas y figuras 4 curvas. Todo escribe dentro de qgrad/. Uso (desde qgrad/):   bash run_17r.sh
set -euo pipefail
cd "$(dirname "$0")"
SA=/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis
R=relax300K/restart.glass300K_N65536_qgrad
[ -f "$R" ] || { echo "Todavia no existe $R -> la relajacion (job 43581066) no termino. Revisar: squeue -u \$USER ; tail relax300K/log.relax.*"; exit 1; }
echo "17r $PWD/$R" > manifest_17r.txt
grep -q '^17r ' manifest_300K.txt || echo "17r $PWD/$R" >> manifest_300K.txt
J1=$(sbatch --parsable --job-name=qg_dump17r batch_qgrad.sbatch dump manifest_17r.txt)
J2=$(sbatch --parsable --job-name=qg_boo17r  batch_qgrad.sbatch boo  manifest_17r.txt)
cd etapa4; J3=$(sbatch --parsable --job-name=qg_rdf17r batch_rdf_qgrad.sbatch ../manifest_17r.txt); cd ..
# rdf: el manifest_rdf.txt ya no se usa; batch_rdf_qgrad.sbatch lee la ruta relativa a etapa4/ -> pasamos ruta absoluta
J4=$(sbatch --parsable --dependency=afterok:$J1 --job-name=qg_voro17r run_voronoi.sbatch 06 07 08 09 10 11 12 13 14 15 16 17u 17r)
cat > laves/run_features_17r.sbatch <<'EOF'
#!/bin/bash
#SBATCH -J qg_feat17r
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge; module load anaconda/2024.02-py311
export PYTHONDONTWRITEBYTECODE=1
SA=/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis
python laves_features.py qgrad 17r $SA/qgrad/snaps/17r.dump $SA/qgrad/boo/17r.boo.dump
EOF
cd laves; J5=$(sbatch --parsable --dependency=afterok:$J1:$J2 run_features_17r.sbatch); cd ..
cat > post_17r.sbatch <<'EOF'
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
EOF
J6=$(sbatch --parsable --dependency=afterok:$J1:$J2:$J3:$J4:$J5 post_17r.sbatch)
echo "encadenados: dump=$J1 boo=$J2 rdf=$J3 voro=$J4 feat=$J5 post=$J6"
echo "Al terminar: repetir checksums -> sbatch checksums.sbatch despues ; diff checksums_antes.txt checksums_despues.txt"
