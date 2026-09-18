#!/usr/bin/env python3
"""Tarea H, pasos 3 y 4 -- merge de step1 (BOO ico) + step2 (Voronoi ico,
Z14/Z15/Z16) en una sola tabla, y validacion cualitativa en el config a
300K de cada tasa contra el patron del piloto (2e10, Tareas F/G):
  - Cu domina ico (BOO y Voronoi)
  - Zr domina Z15 y Z16
  - Ag intermedio en Z14 (entre Zr y Cu)

Si alguna tasa rompe el patron, el script PARA (assert) antes de escribir
la tabla final -- no la da por buena silenciosamente.
"""
import csv

RATES = ["2e10", "1e11", "1e12"]
SPECIES = ["Zr", "Cu", "Ag"]


def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


boo = load_csv("step1_boo_ico.csv")
voro = load_csv("step2_voronoi.csv")

boo_idx = {(r["tasa"], r["config"], r["especie"]): r for r in boo}
voro_idx = {(r["tasa"], r["config"], r["especie"]): r for r in voro}

assert set(boo_idx) == set(voro_idx), "step1 y step2 no tienen las mismas claves"

merged = []
for key in sorted(boo_idx):
    b = boo_idx[key]
    v = voro_idx[key]
    merged.append({
        "tasa": b["tasa"], "config": b["config"], "T_nominal": b["T_nominal"],
        "especie": b["especie"], "N_total_especie": b["N_total_especie"],
        "N_ico_BOO": b["N_ico_BOO"], "N_ico_Voronoi": v["N_ico_Voronoi"],
        "N_Z14": v["N_Z14"], "N_Z15": v["N_Z15"], "N_Z16": v["N_Z16"],
    })

# --------------------------------------------------------------- paso 3
print("=== paso 3: validacion cualitativa a 300K (config 17) por tasa ===")
problems = []
for rate in RATES:
    row_by_sp = {r["especie"]: r for r in merged
                 if r["tasa"] == rate and r["config"] == "17"}
    ico_boo = {sp: int(row_by_sp[sp]["N_ico_BOO"]) for sp in SPECIES}
    ico_vor = {sp: int(row_by_sp[sp]["N_ico_Voronoi"]) for sp in SPECIES}
    z14 = {sp: int(row_by_sp[sp]["N_Z14"]) for sp in SPECIES}
    z15 = {sp: int(row_by_sp[sp]["N_Z15"]) for sp in SPECIES}
    z16 = {sp: int(row_by_sp[sp]["N_Z16"]) for sp in SPECIES}

    print(f"\n[{rate}] T=300K")
    print(f"  ico BOO:     Zr={ico_boo['Zr']:5d}  Cu={ico_boo['Cu']:5d}  Ag={ico_boo['Ag']:5d}")
    print(f"  ico Voronoi: Zr={ico_vor['Zr']:5d}  Cu={ico_vor['Cu']:5d}  Ag={ico_vor['Ag']:5d}")
    print(f"  Z14:         Zr={z14['Zr']:5d}  Cu={z14['Cu']:5d}  Ag={z14['Ag']:5d}")
    print(f"  Z15:         Zr={z15['Zr']:5d}  Cu={z15['Cu']:5d}  Ag={z15['Ag']:5d}")
    print(f"  Z16:         Zr={z16['Zr']:5d}  Cu={z16['Cu']:5d}  Ag={z16['Ag']:5d}")

    checks = {
        "Cu domina ico (BOO)": ico_boo["Cu"] == max(ico_boo.values()),
        "Cu domina ico (Voronoi)": ico_vor["Cu"] == max(ico_vor.values()),
        "Zr domina Z15": z15["Zr"] == max(z15.values()),
        "Zr domina Z16": z16["Zr"] == max(z16.values()),
        "Ag intermedio en Z14 (entre Zr y Cu)":
            min(z14["Zr"], z14["Cu"]) <= z14["Ag"] <= max(z14["Zr"], z14["Cu"]),
    }
    for label, ok in checks.items():
        mark = "OK" if ok else "FALLA"
        print(f"    [{mark}] {label}")
        if not ok:
            problems.append((rate, label))

if problems:
    print("\n" + "=" * 60)
    print("PROBLEMAS DETECTADOS -- no se escribe la tabla final:")
    for rate, label in problems:
        print(f"  [{rate}] {label}")
    raise SystemExit(1)

print("\nTodas las tasas mantienen el patron cualitativo del piloto (2e10).")

# --------------------------------------------------------------- paso 4
import os
os.makedirs(".", exist_ok=True)
out = "tabla_clusters_51.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(merged[0].keys()))
    w.writeheader()
    w.writerows(merged)
print(f"\n-> {out} ({len(merged)} filas = {len(RATES)} tasas x 17 configs x "
      f"{len(SPECIES)} especies)")
