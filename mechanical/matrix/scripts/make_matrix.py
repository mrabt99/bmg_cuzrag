#!/usr/bin/env python3
"""Crea matrix/<sample>_<T>K_<erate>_<dir>/ (copia de in.tension + sbatch) para todas las etapas y escribe/actualiza manifest.csv.
No lanza nada. Idempotente: no pisa casos existentes ni sus estados. Los casos 2e10 y 1e12 en x a 300 K y 1e9 (los pilotos) se ENLAZAN (symlink), no se repiten.
Etapas:  1: 300 K; 1e9 y 1e10; 2e10 1e11 1e12; x y z      2: 450 y 540 K; 1e9; tres tasas; x y z
         3: 300 K; 1e8; tres tasas; x y z (3 h)           4: 750 y 900 K; 1e9; solo 2e10; x y z (por encima de Tg(2e10) = 641.7 K)
         5: 300 K; 1e9, 1e10 y 1e8; qgrad y qgrad_ext (perfiles graduados, sin tasa unica); x y z
         extra: Cij a 300 K y a 0 K de 1e11, 1e12, qgrad y qgrad_ext (mismo procedimiento que 2e10, copiando mechanical/elastic_2e10)
Uso: python make_matrix.py"""
import os, csv, shutil, re
M = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))          # .../mechanical/matrix
MECH = os.path.dirname(M)
SAMPLES = ["2e10", "1e11", "1e12"]; DIRS = ["x", "y", "z"]
STAGES = {1: dict(T=[300], erate=["1e9", "1e10"], samples=SAMPLES), 2: dict(T=[450, 540], erate=["1e9"], samples=SAMPLES),
          3: dict(T=[300], erate=["1e8"], samples=SAMPLES), 4: dict(T=[750, 900], erate=["1e9"], samples=["2e10"]),
          5: dict(T=[300], erate=["1e9", "1e10", "1e8"], samples=["qgrad", "qgrad_ext"])}
LINKED = {("2e10", 300, "1e9", "x"): "../pilot_tension", ("1e12", 300, "1e9", "x"): "../pilot_tension_1e12"}
NHOLD = 50000; DT_PS = 0.001
def steps(T, erate):
    nramp = round((T - 300) / 1e12 / 1e-15); nsteps = round(0.2 / (float(erate) * 1e-12 * DT_PS)); return nramp + NHOLD + nsteps
def walltime(n):          # a ~300 pasos/s medidos con 48 ranks (piloto: 250000 pasos = 14 min); margen >= 2x salvo etapa 3 (3 h pedidas, estimado 114 min)
    return "00:30:00" if n <= 100_000 else "01:00:00" if n <= 300_000 else "01:30:00" if n <= 600_000 else "02:00:00" if n <= 900_000 else "03:00:00"
SBATCH = """#!/bin/bash
# Traccion: muestra {sample}, T = {T} K, erate = {erate} s^-1, eje {dir}. Pasos estimados: {steps} (~{minutes} min a 300 pasos/s). standby (<4 h), pmi2, 48 ranks.
#SBATCH --job-name=t_{sample}_{T}_{erate}_{dir}
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err
#SBATCH -A strachan
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time={walltime}
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load gcc/12.2.0 openmpi/4.1.4 lammps/20220623
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in in.tension -var sample {sample} -var T {T} -var erate {erate} -var dir {dir} \\
     -log log.tens_{sample}_{T}K_{erate}_{dir}.${{SLURM_JOB_ID}}
"""
old = {}
mpath = os.path.join(M, "manifest.csv")
if os.path.exists(mpath): old = {r["case_id"]: r for r in csv.DictReader(open(mpath))}
rows = []
def add(r):
    if r["case_id"] in old: r["status"], r["jobid"] = old[r["case_id"]]["status"], old[r["case_id"]]["jobid"]
    rows.append(r)
for st, g in STAGES.items():
    for erate in g["erate"]:
        for T in g["T"]:
            for sample in g["samples"]:
                for d in DIRS:
                    cid = f"{sample}_{T}K_{erate}_{d}"; path = cid; n = steps(T, erate); wt = walltime(n)
                    r = dict(case_id=cid, stage=st, sample=sample, T=T, erate=erate, dir=d, path=path, status="pendiente", jobid="", steps_est=n, walltime=wt, note="")
                    if (sample, T, erate, d) in LINKED:
                        tgt = LINKED[(sample, T, erate, d)]; r.update(path=cid, status="ok (enlazado)", walltime="", note=f"piloto ya corrido: symlink a {tgt[3:]}")
                        if not os.path.lexists(os.path.join(M, cid)): os.symlink(tgt, os.path.join(M, cid))
                    else:
                        p = os.path.join(M, cid)
                        if "--refresh-inputs" in __import__("sys").argv and os.path.isdir(p) and old.get(cid, {}).get("status", "pendiente") == "pendiente":
                            shutil.copy(os.path.join(M, "in.tension"), os.path.join(p, "in.tension"))     # re-copia el input maestro SOLO a casos aun pendientes
                        if not os.path.exists(p):
                            os.makedirs(p); shutil.copy(os.path.join(M, "in.tension"), os.path.join(p, "in.tension"))
                            open(os.path.join(p, f"slurm_tension.sbatch"), "w").write(SBATCH.format(sample=sample, T=T, erate=erate, dir=d, steps=n, minutes=round(n / 300 / 60), walltime=wt))
                    add(r)
# ---- extra: elasticas de 1e11 y 1e12 (copiando elastic_2e10 y cambiando el restart de entrada)
SRC = os.path.join(MECH, "elastic_2e10"); RST = {"1e11": "rate_1e11/prod_1e11/restart.glass300K_N65536_1e11", "1e12": "rate_1e12_reference/restart.glass300K_N65536_1e12",
       "qgrad": "structural_analysis/qgrad/relax300K/restart.glass300K_N65536_qgrad", "qgrad_ext": "qgrad_ext/restart.glass300K_N65536_gradext"}
FILES = {"300K": ["init.mod", "potential.mod", "in.elastic_prep", "in.elastic_dir", "run_all.sh", "slurm_elastic300K.sbatch"],
         "0K": ["init.mod", "potential.mod", "in.elastic_prep", "in.elastic_dir", "run_0K.sh", "slurm_elastic0K.sbatch"]}
for sample in ("1e11", "1e12", "qgrad", "qgrad_ext"):
    for tt in ("300K", "0K"):
        dst = os.path.join(MECH, f"elastic_{sample}", tt)
        if not os.path.exists(dst):
            os.makedirs(dst)
            for f in FILES[tt]:
                t = open(os.path.join(SRC, tt, f)).read()
                if f == "init.mod":
                    t, k = re.subn(r"(variable rst\s+index )\S+", r"\g<1>/scratch/negishi/aaabt/bmg_cuzrag/" + RST[sample], t); assert k == 1
                if f.startswith("slurm_"): t = t.replace("--job-name=elastic300K", f"--job-name=el300_{sample}").replace("--job-name=elastic0K", f"--job-name=el0K_{sample}")
                open(os.path.join(dst, f), "w").write(t); os.chmod(os.path.join(dst, f), os.stat(os.path.join(SRC, tt, f)).st_mode)
        cid = f"elastic_{sample}_{tt}"
        add(dict(case_id=cid, stage="extra", sample=sample, T=tt.replace("K", ""), erate="-", dir="-", path=f"../elastic_{sample}/{tt}", status="pendiente", jobid="",
                 steps_est="", walltime="03:00:00" if tt == "300K" else "01:30:00", note="Cij (mismo procedimiento que elastic_2e10)"))
hdr = ["case_id", "stage", "sample", "T", "erate", "dir", "path", "status", "jobid", "steps_est", "walltime", "note"]
with open(mpath, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=hdr); w.writeheader(); w.writerows(rows)
from collections import Counter
print(f"manifest.csv: {len(rows)} casos"); print({k: v for k, v in sorted(Counter(str(r['stage']) for r in rows).items())})
