#!/usr/bin/env python3
"""Actualiza manifest.csv: para los casos 'corriendo' consulta sacct y el log. Traccion: ok si COMPLETED y el log tiene '### TRACCION TERMINADA'.
Elasticas: ok si COMPLETED y existe la ultima columna (Ccol6.txt, o up_5e-4/Ccol6.txt a 0 K). Cualquier otro estado final -> 'fallo'. Si la traccion termino pero check_integrity.py encuentra NUL/cuadros faltantes -> 'danado'. Uso: python update_status.py"""
import os, csv, glob, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_integrity import check
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); mpath = os.path.join(M, "manifest.csv")
rows = list(csv.DictReader(open(mpath))); ch = 0
for r in rows:
    if r["status"] != "corriendo" or not r["jobid"]: continue
    st = subprocess.run(["sacct", "-j", r["jobid"], "-X", "-n", "-P", "-o", "State"], capture_output=True, text=True).stdout.strip().split("\n")[0].split()[0] if True else ""
    p = os.path.join(M, r["path"])
    if st in ("PENDING", "RUNNING", "COMPLETING", "REQUEUED", ""): continue
    if st == "COMPLETED":
        if r["stage"] == "extra": good = os.path.exists(os.path.join(p, "up_5e-4/Ccol6.txt" if r["T"] == "0" else "Ccol6.txt"))
        else: good = any("### TRACCION TERMINADA" in open(l, errors="ignore").read() for l in glob.glob(os.path.join(p, "log.tens_*")))
        r["status"] = "ok" if good else "fallo"
        if good and r["stage"] != "extra":
            pr = check(p)
            if pr: r["status"] = "danado"; r["note"] = (r["note"] + " " if r["note"] else "") + "INTEGRIDAD: " + "; ".join(pr)
    else: r["status"] = "fallo"
    r["note"] = (r["note"] + " " if r["note"] else "") + f"sacct={st}"; ch += 1; print(f"{r['case_id']}: {r['status']} ({st})")
with open(mpath, "w", newline="") as f: w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"{ch} cambios;", {k: sum(1 for r in rows if r['status'] == k) for k in sorted({r['status'] for r in rows})})
