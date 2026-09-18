#!/usr/bin/env python3
"""Tarea H, paso 1 -- criterio de ico por BOO sobre los 51 *.boo.dump ya
calculados (structural_analysis/boo/<tasa>/NN.boo.dump), SOLO LECTURA,
no se toca boo/.

Config NN (01..17) -> T nominal: mismo orden de manifest.txt en las 3
tasas (1900K bajando de 100 en 100 hasta 300K relajado), confirmado al
construir los manifiestos en la Tarea E/G: NN=1->1900K, ..., NN=17->300K.
"""
import csv
import sys

sys.path.insert(0, "../classify")
from boo_types import load_boo_dump, mask_ico, SPECIES  # noqa: E402

RATES = ["2e10", "1e11", "1e12"]
NCONFIGS = 17


def t_nominal(nn):
    return 1900 - (nn - 1) * 100


rows = []
for rate in RATES:
    for nn in range(1, NCONFIGS + 1):
        path = f"../boo/{rate}/{nn:02d}.boo.dump"
        data = load_boo_dump(path)
        types = data[:, 1].astype(int)
        ico = mask_ico(data)
        T = t_nominal(nn)
        for t, name in SPECIES.items():
            m = types == t
            rows.append({
                "tasa": rate, "config": f"{nn:02d}", "T_nominal": T,
                "especie": name,
                "N_total_especie": int(m.sum()),
                "N_ico_BOO": int((ico & m).sum()),
            })
        print(f"[{rate}] {nn:02d} (T={T}K): ico BOO total="
              f"{ico.sum()} -- " + ", ".join(
                  f"{name}={(ico & (types == t)).sum()}"
                  for t, name in SPECIES.items()))

with open("step1_boo_ico.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"\n-> step1_boo_ico.csv ({len(rows)} filas = {len(RATES)} tasas x "
      f"{NCONFIGS} configs x {len(SPECIES)} especies)")
