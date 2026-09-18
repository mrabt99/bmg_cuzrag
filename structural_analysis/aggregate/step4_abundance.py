#!/usr/bin/env python3
"""Tarea I - aggregate_abundance. Puramente sobre tabla_clusters_51.csv,
no recalcula nada de BOO/Voronoi.

Por tasa, por especie, por T: fraccion = N_tipo / N_total_especie, para
n_ico_boo, n_ico_voronoi, n_Z14, n_Z15, n_Z16.

n_ico_BOO y n_ico_Voronoi se guardan como columnas separadas a proposito
(decision del usuario): difieren en magnitud (~4-5x) y en sensibilidad a
la tasa, mezclarlos ocultaria eso.
"""
import csv

RATES = ["2e10", "1e11", "1e12"]

with open("tabla_clusters_51.csv") as f:
    rows = list(csv.DictReader(f))

for rate in RATES:
    out_rows = []
    for r in rows:
        if r["tasa"] != rate:
            continue
        ntot = int(r["N_total_especie"])
        out_rows.append({
            "T": int(r["T_nominal"]),
            "especie": r["especie"],
            "n_ico_boo": int(r["N_ico_BOO"]) / ntot,
            "n_ico_voronoi": int(r["N_ico_Voronoi"]) / ntot,
            "n_Z14": int(r["N_Z14"]) / ntot,
            "n_Z15": int(r["N_Z15"]) / ntot,
            "n_Z16": int(r["N_Z16"]) / ntot,
        })
    out_rows.sort(key=lambda x: (-x["T"], x["especie"]))
    out_path = f"abundancia_{rate}.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["T", "especie", "n_ico_boo",
                                          "n_ico_voronoi", "n_Z14", "n_Z15", "n_Z16"])
        w.writeheader()
        for row in out_rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v)
                       for k, v in row.items()})
    print(f"-> {out_path} ({len(out_rows)} filas)")
