#!/usr/bin/env python3
"""Relanza los casos 'danado' o 'fallo' del manifest con un TOPE de jobs simultaneos (por defecto 6 = 288 cores): reparte los casos en N carreras (cadenas) y
cada job depende (afterany) del anterior de su carrera, asi nunca corren mas de N a la vez aunque el planificador lo permitiria. Las salidas del intento anterior
se MUEVEN a <caso>/intento<k>_<jobid>/ (no se borran). Por defecto ensayo en seco; con --go envia. Uso: python relaunch.py [--max-concurrent 6] [--go]"""
import os, sys, csv, glob, shutil, subprocess
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); mpath = os.path.join(M, "manifest.csv")
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
N = int(arg("--max-concurrent", 6)); GO = "--go" in sys.argv
rows = list(csv.DictReader(open(mpath))); todo = [r for r in rows if r["status"] in ("danado", "fallo")]
def hours(r):
    h, m, s = r["walltime"].split(":"); return int(h) + int(m) / 60
KEEP = ("in.tension", "slurm_tension.sbatch")
lanes = [[] for _ in range(N)]; load = [0.0] * N
for r in sorted(todo, key=lambda r: -hours(r)):
    i = load.index(min(load)); lanes[i].append(r); load[i] += hours(r)
print(f"{len(todo)} casos a relanzar en {N} carreras (tope {N} jobs = {48*N} cores); carga por carrera (h de walltime): {[round(x,1) for x in load]}" + ("" if GO else "   [ENSAYO EN SECO]"))
for i, lane in enumerate(lanes):
    print(f"  carrera {i+1}: " + " -> ".join(f"{r['case_id']}({r['walltime'][:4]})" for r in lane))
if GO:
    for lane in lanes:
        prev = None
        for r in lane:
            p = os.path.join(M, r["path"])
            if r["stage"] != "extra":                       # mover salidas viejas a intentoK_<jobid>/
                old = [f for f in os.listdir(p) if f not in KEEP and not os.path.isdir(os.path.join(p, f))]
                if old:
                    k = 1 + len(glob.glob(os.path.join(p, "intento*_*"))); d = os.path.join(p, f"intento{k}_{r['jobid']}"); os.makedirs(d)
                    for f in old: shutil.move(os.path.join(p, f), os.path.join(d, f))
            else:
                old = [f for f in os.listdir(p) if f not in ("init.mod", "potential.mod", "in.elastic_prep", "in.elastic_dir", "run_all.sh", "run_0K.sh") and not f.startswith("slurm_elastic") or f.endswith((".out", ".err"))]
                old = [f for f in old if os.path.isfile(os.path.join(p, f)) and not f.endswith(".sbatch")]
                if old:
                    d = os.path.join(p, f"intento1_{r['jobid']}"); os.makedirs(d)
                    for f in old: shutil.move(os.path.join(p, f), os.path.join(d, f))
            sb = next(f for f in os.listdir(p) if f.startswith("slurm_") and f.endswith(".sbatch"))
            cmd = ["sbatch", "--parsable"] + ([f"--dependency=afterany:{prev}"] if prev else []) + [sb]
            out = subprocess.run(cmd, cwd=p, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            if out.returncode: print("  ERROR", r["case_id"], out.stderr.strip()); break
            jid = out.stdout.strip().split(";")[0]; r["note"] = (r["note"] + " | " if r["note"] else "") + f"relanzado (antes job {r['jobid']}, {r['status']})"; r["jobid"] = jid; r["status"] = "corriendo"; prev = jid
            print(f"  {r['case_id']:26s} job {jid}" + (f" (despues de {out.args[2].split(':')[1]})" if prev and len(out.args) > 3 and out.args[2].startswith('--dependency') else ""))
    with open(mpath, "w", newline="") as f: w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
