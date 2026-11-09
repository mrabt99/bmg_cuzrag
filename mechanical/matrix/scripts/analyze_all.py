#!/usr/bin/env python3
"""Analiza cada caso de traccion 'ok' (curva con analyze_curve.py; localizacion con localization.py via OVITO/anaconda) y arma matrix/resumen.csv:
case_id, stage, sample, T, erate, dir, E, sigma_max, eps_sigma_max, sigma_flujo, delta_sigma, psi_max, psi_w00_02 ... psi_w18_20 (psi por ventana de 0.02).
Uso (con el python del venv_structural): python analyze_all.py [--force] [--stage N] [--only case_id]   Los casos ya analizados (y no mas viejos que su curva) se saltean."""
import os, sys, csv, glob, json, subprocess
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); S = os.path.dirname(os.path.abspath(__file__))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FORCE = "--force" in sys.argv; only = arg("--only"); stage = arg("--stage")
rows = [r for r in csv.DictReader(open(os.path.join(M, "manifest.csv"))) if r["status"].startswith("ok") and r["stage"] != "extra"]
if only: rows = [r for r in rows if r["case_id"] == only]
if stage: rows = [r for r in rows if r["stage"] == stage]
PYV = sys.executable
def run_loc(p, d, sample):
    sample = sample if sample in ("2e10", "1e11", "1e12") else "2e10"   # cutoff comun 3.75 A (primer minimo de g(r)); qgrad* no tiene g(r) en rdf_all.csv
    cmd = f"module purge; module load anaconda/2024.02-py311 >/dev/null 2>&1; export PYTHONDONTWRITEBYTECODE=1; python {S}/localization.py --dir {p} --axis {d} --rate {sample} --out {p}/analysis/localization"
    return subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
out = []
for r in rows:
    p = os.path.realpath(os.path.join(M, r["path"])); curve = glob.glob(f"{p}/*.curve.dat")[0]
    cj, pc = f"{p}/analysis/curva_resumen.json", f"{p}/analysis/localization/psi_ventanas.csv"
    fresh = all(os.path.exists(f) and os.path.getmtime(f) >= os.path.getmtime(curve) for f in (cj, pc))
    if FORCE or not fresh:
        a = subprocess.run([PYV, f"{S}/analyze_curve.py", curve, "--dir", r["dir"], "--label", r["case_id"]], capture_output=True, text=True)
        if a.returncode: print(r["case_id"], "analyze_curve ERROR:", a.stderr[-300:]); continue
        l = run_loc(p, r["dir"], r["sample"])
        if not os.path.exists(pc): print(r["case_id"], "localization ERROR:", (l.stderr or l.stdout)[-400:]); continue
        print("analizado:", r["case_id"], flush=True)
    c = json.load(open(cj)); ps = [float(x["psi_std_eta"]) for x in csv.DictReader(open(pc))]
    d = dict(case_id=r["case_id"], stage=r["stage"], sample=r["sample"], T=r["T"], erate=r["erate"], dir=r["dir"], E_GPa=f"{c['E']:.2f}", sigma_max_GPa=f"{c['sigma_max']:.4f}",
             eps_sigma_max=f"{c['eps_sigma_max']:.4f}", sigma_flujo_GPa=f"{c['sigma_flow']:.4f}", delta_sigma_GPa=f"{c['delta_sigma']:.4f}", psi_max=f"{max(ps):.5f}")
    for k, v in enumerate(ps): d[f"psi_w{2*k:02d}_{2*k+2:02d}"] = f"{v:.5f}"
    out.append(d)
if out:
    mp = os.path.join(M, "resumen.csv"); prev = {x["case_id"]: x for x in csv.DictReader(open(mp))} if os.path.exists(mp) else {}
    for d in out: prev[d["case_id"]] = d
    hdr = list(out[0].keys()); allr = sorted(prev.values(), key=lambda x: (int(x["stage"]), x["sample"], float(x["T"]), x["erate"], x["dir"]))
    with open(mp, "w", newline="") as f: w = csv.DictWriter(f, fieldnames=hdr); w.writeheader(); w.writerows(allr)
    print(f"resumen.csv: {len(allr)} casos")
