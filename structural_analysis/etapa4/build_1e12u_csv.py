#!/usr/bin/env python3
"""Arma los CSV de RDF de la 1e12 UNIFICADA (rate = 1e12u) con el MISMO formato que los de las otras tasas,
pero en archivos aparte (sufijo _1e12u); los CSV existentes no se tocan (17/11/2026).

Entrada : rdf/raw/1e12u, rdf_partial/raw/1e12u, rdf_sametype/raw/1e12u, rdf_cross/raw/1e12u (NN.rdf.txt + NN.log)
Salida  : rdf/rdf_all_1e12u.csv, rdf/rdf_qc_1e12u.csv, rdf_partial/rdf_partial_all_1e12u.csv,
          rdf_sametype/rdf_sametype_all_1e12u.csv, rdf_cross/rdf_cross_all_1e12u.csv
Chequeos: g_total de las 4 variantes idéntico; g_ZrCu == g_CuZr; N(<9,99 Å) / ideal ≈ 1.
Uso (desde etapa4/): PYTHONDONTWRITEBYTECODE=1 /scratch/negishi/aaabt/venv_structural/bin/python build_1e12u_csv.py
"""
import csv, re
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RATE, NBIN = "1e12u", 500
VARIANTES = {
    "rdf": (["g_total", "Ncum"], "rdf_all"),
    "rdf_partial": (["g_total", "Ncum_total", "g_Zr", "Ncum_Zr", "g_Cu", "Ncum_Cu"], "rdf_partial_all"),
    "rdf_sametype": (["g_total", "Ncum_total", "g_ZrZr", "Ncum_ZrZr", "g_CuCu", "Ncum_CuCu"], "rdf_sametype_all"),
    "rdf_cross": (["g_total", "Ncum_total", "g_ZrCu", "Ncum_ZrCu", "g_CuZr", "Ncum_CuZr"], "rdf_cross_all"),
}


def T_of(nn):
    return 1900 - 100 * (nn - 1)


def leer(path, ncol):
    a = np.array([[float(x) for x in l.split()[1:]] for l in open(path)
                  if not l.startswith("#") and len(l.split()) == ncol + 2])
    assert a.shape == (NBIN, ncol + 1), (path, a.shape)
    return a


gtot = {}
qc = []
for var, (cols, base) in VARIANTES.items():
    filas = []
    for nn in range(1, 18):
        a = leer(HERE / var / "raw" / RATE / f"{nn:02d}.rdf.txt", len(cols))
        T = T_of(nn)
        gtot.setdefault(T, []).append(a[:, 1])
        if var == "rdf_cross":
            assert np.array_equal(a[:, 3], a[:, 5]), f"g_ZrCu != g_CuZr en {nn}"
        if var == "rdf":
            m = re.search(r"QC_BOX lx=(\S+) ly=\S+ lz=\S+ natoms=(\d+) numdens=(\S+)",
                          (HERE / var / "raw" / RATE / f"{nn:02d}.log").read_text())
            L, rho = float(m.group(1)), float(m.group(3))
            assert int(m.group(2)) == 65536
            ideal = rho * 4 / 3 * np.pi * a[-1, 0] ** 3
            qc.append([RATE, T, f"{nn:02d}", L, rho, a[-1, 2], ideal, a[-1, 2] / ideal])
        for row in a:
            filas.append([RATE, T, f"{row[0]:.4f}"] + [f"{x:.6f}" for x in row[1:]])
    with open(HERE / var / f"{base}_{RATE}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["rate", "T", "r"] + cols); w.writerows(filas)
    print(f"{var}: {len(filas)} filas -> {var}/{base}_{RATE}.csv")
with open(HERE / "rdf" / f"rdf_qc_{RATE}.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "T", "config", "L", "numdens", "Ncum_rmax", "Ncum_ideal", "ratio"])
    w.writerows([[r[0], r[1], r[2], f"{r[3]:.5f}", f"{r[4]:.6f}", f"{r[5]:.3f}", f"{r[6]:.3f}", f"{r[7]:.4f}"] for r in qc])
dmax = max(np.max(np.abs(np.array(v) - v[0])) for v in gtot.values())
print(f"chequeo g_total igual en las 4 variantes: max|dif| = {dmax:.2e}")
print("chequeo N(<9,99)/ideal: " + f"{min(r[7] for r in qc):.4f} - {max(r[7] for r in qc):.4f}")
