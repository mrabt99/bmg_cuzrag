#!/usr/bin/env python3
"""Lanza una etapa leyendo manifest.csv. Por defecto es un ENSAYO EN SECO (lista lo que enviaria). Solo envia con --go.
Antes de enviar verifica que las etapas anteriores esten ok (salvo --force). Casos 'pendiente' -> sbatch en su carpeta -> status 'corriendo' + jobid.
Uso: python launch_stage.py --stage 1 [--go] [--force] [--max-concurrent N]      (stage: 1 2 3 4 5 extra; --max-concurrent N = a lo sumo N jobs a la vez)"""
import os, sys, csv, subprocess
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); mpath = os.path.join(M, "manifest.csv")
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
stage = arg("--stage"); GO = "--go" in sys.argv; FORCE = "--force" in sys.argv
assert stage in ("1", "2", "3", "4", "5", "extra"), "--stage 1|2|3|4|5|extra"
MAXC = int(sys.argv[sys.argv.index("--max-concurrent") + 1]) if "--max-concurrent" in sys.argv else 0   # tope de jobs simultaneos (0 = sin tope)
rows = list(csv.DictReader(open(mpath))); order = ["1", "2", "3", "4", "5"]
def okrow(r): return r["status"].startswith("ok")
if stage in order and not FORCE:
    prev = [r for r in rows if r["stage"] in order[:order.index(stage)] and not okrow(r)]
    if prev: sys.exit(f"ABORTO: {len(prev)} casos de etapas anteriores no estan ok (p. ej. {prev[0]['case_id']}: {prev[0]['status']}). Use --force para ignorar.")
todo = [r for r in rows if r["stage"] == stage and r["status"] == "pendiente"]
print(f"etapa {stage}: {len(todo)} casos pendientes" + ("" if GO else "  [ENSAYO EN SECO: no se envia nada; agregar --go]"))
def hrs(r): h, m, s = r["walltime"].split(":"); return int(h) + int(m) / 60
lanes = [[] for _ in range(MAXC or len(todo) or 1)]; load = [0.0] * len(lanes)       # reparto LPT en carreras; con --max-concurrent cada job depende (afterany) del anterior de su carrera
for r in sorted(todo, key=lambda r: -hrs(r)):
    i = load.index(min(load)); lanes[i].append(r); load[i] += hrs(r)
if MAXC: print(f"  tope de {MAXC} jobs simultaneos ({48*MAXC} cores): {len([l for l in lanes if l])} carreras, carga (h de walltime) {[round(x, 1) for x in load]}")
for lane in lanes:
    prev = None
    for r in lane:
        p = os.path.join(M, r["path"]); sb = next(f for f in os.listdir(p) if f.startswith("slurm_") and f.endswith(".sbatch"))
        if not GO: print(f"  {r['case_id']:28s} wall {r['walltime']}  ({sb} en {r['path']})"); continue
        cmd = ["sbatch", "--parsable"] + ([f"--dependency=afterany:{prev}"] if prev else []) + [sb]
        out = subprocess.run(cmd, cwd=p, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if out.returncode: print(f"  {r['case_id']}: ERROR sbatch: {out.stderr.strip()}"); break
        r["jobid"] = out.stdout.strip().split(";")[0]; r["status"] = "corriendo"; prev = r["jobid"]; print(f"  {r['case_id']:28s} job {r['jobid']}" + (f"  (tras {cmd[2].split(':')[1]})" if len(cmd) > 3 else ""))
if GO:
    with open(mpath, "w", newline="") as f: w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
