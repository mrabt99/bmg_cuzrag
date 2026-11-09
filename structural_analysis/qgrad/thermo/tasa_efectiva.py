#!/usr/bin/env python3
"""HIPOTESIS (no resultado): 'tasa efectiva' de qgrad segun Ep(300 K) y lo que predeciria para los clusters.
SUPUESTO CLAVE: Ep(300 K vidrio relajado) es LINEAL en log10(tasa) (se verifica con 3 tasas, ver pendientes) y la abundancia de cada cluster
es log-lineal en log10(tasa) (ln N lineal). Ep = run 0 de restart.glass300K_* (ep_run0_resultados.txt); pools = tabla_regimen_qgrad.csv (vitreo 500/400/300 K).
Salida: tasa_efectiva.csv"""
import csv, math, re
import numpy as np
EP = {m.group(1): float(m.group(2)) for l in open("ep_run0_resultados.txt") for m in [re.search(r"tag=(\S+) .*Ep_at=(\S+)", l)] if m}
RATE = {"2e10": 2e10, "1e11": 1e11, "1e12": 1e12}; LG = {k: math.log10(v) for k, v in RATE.items()}
print("Ep run 0 (eV/at):", {k: round(v, 6) for k, v in EP.items()})
print("\nPendientes dEp/dlog10(tasa) (meV/at por decada):")
sl = {"2e10->1e11": 1000 * (EP["1e11"] - EP["2e10"]) / (LG["1e11"] - LG["2e10"]), "1e11->1e12": 1000 * (EP["1e12"] - EP["1e11"]) / (LG["1e12"] - LG["1e11"]),
      "2e10->1e12": 1000 * (EP["1e12"] - EP["2e10"]) / (LG["1e12"] - LG["2e10"])}
for k, v in sl.items(): print(f"  {k}: {v:.3f}")
dEp = 1000 * (EP["qgrad"] - EP["2e10"]); print(f"\nEp(qgrad) - Ep(2e10) = {dEp:+.3f} meV/at")
dlog = {k: dEp / v for k, v in sl.items()}
print("Delta log10(tasa) efectiva respecto de 2e10 (decadas):", {k: round(v, 3) for k, v in dlog.items()})
print(f"  -> tasa efectiva = 2e10 * 10^{dlog['2e10->1e11']:.3f} = {2e10 * 10 ** dlog['2e10->1e11']:.3e} K/s  (pendiente 2e10->1e11; con 1e11->1e12: {2e10 * 10 ** dlog['1e11->1e12']:.3e})")
# pools vitreos
P = {}
for r in csv.DictReader(open("../aggregate/tabla_regimen_qgrad.csv")):
    if r["regimen"] == "vitreo" and r["N_pool_qgrad"]:
        P[(r["referencia"], r["especie"], r["tipo_cluster"])] = (int(r["N_pool_qgrad"]), int(r["N_pool_ref"]))
rows = []
print(f"\n{'cluster':16s} {'medido qgrad/2e10':>20s} {'pred (2 pts 2e10-1e11)':>24s} {'pred (ajuste 3 tasas)':>22s} {'medido-pred (2 pts)':>20s}")
for sp, tc, lab in [("Cu", "ico_boo", "Cu ico BOO"), ("Cu", "ico_voronoi", "Cu ico Voronoi"), ("Zr", "Z14_voronoi", "Zr Z14"), ("Zr", "Z15_voronoi", "Zr Z15"), ("Zr", "Z16_voronoi", "Zr Z16")]:
    nq, n21 = P[("2e10", sp, tc)]; _, n11 = P[("1e11", sp, tc)]; _, n12 = P[("1e12", sp, tc)]
    meas = nq / n21; e_meas = meas * math.sqrt(1 / nq + 1 / n21)
    ex = -dlog["2e10->1e11"] / (LG["1e11"] - LG["2e10"])            # >0: qgrad esta 'ex' de camino (en unidades del salto 2e10->1e11) hacia mas lento
    r21 = n21 / n11; pred2 = r21 ** ex                                # prediccion con 2 puntos (2e10 y 1e11), ln N lineal en log10(tasa)
    e_pred2 = pred2 * abs(ex) * math.sqrt(1 / n21 + 1 / n11)          # solo Poisson de los dos pools de referencia (piso)
    x = np.array([LG["1e12"], LG["1e11"], LG["2e10"]]); y = np.log([n12, n11, n21]); bslope, _ = np.polyfit(x, y, 1)
    pred3 = math.exp(bslope * dlog["2e10->1e11"])                     # idem con ajuste lineal de ln N vs log10(tasa) a las 3 tasas
    z = (meas - pred2) / math.hypot(e_meas, e_pred2)
    print(f"{lab:16s} {meas:12.3f} +- {e_meas:.3f} {pred2:16.3f} +- {e_pred2:.3f} {pred3:18.3f} {meas - pred2:+14.3f} ({z:+.1f}σ)")
    rows.append([lab, f"{meas:.4f}", f"{e_meas:.4f}", f"{pred2:.4f}", f"{e_pred2:.4f}", f"{pred3:.4f}", f"{z:.2f}"])
with open("tasa_efectiva.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["cluster", "razon_medida_qgrad/2e10", "error", "razon_predicha_2pts", "error_pred", "razon_predicha_ajuste3", "(medida-pred)/error"]); w.writerows(rows)
print("\nSensibilidad de la prediccion Cu-ico BOO a Delta log10(tasa) = 0.10 / 0.15 / 0.20 decadas (2 pts):")
nq, n21 = P[("2e10", "Cu", "ico_boo")]; _, n11 = P[("1e11", "Cu", "ico_boo")]
for d in (0.10, 0.15, 0.20): print(f"   {d:.2f} dec -> razon predicha {(n21 / n11) ** (d / (LG['1e11'] - LG['2e10'])):.3f}")
# como evoluciona la razon qgrad/2e10 con T (antes y despues del tramo rapido)
N = {}
for r in csv.DictReader(open("../aggregate/tabla_compare_rates_4curvas.csv")):
    if r["especie"] == "Cu" and r["tipo_cluster"] in ("ico_boo", "ico_voronoi") or (r["especie"] == "Zr" and r["tipo_cluster"] == "Z16_voronoi"):
        c = r["config_qgrad"]
        if r["rate"] in ("2e10",) or (r["rate"] == "qgrad" and c != "17u"): N[(r["rate"], int(r["T"]), r["especie"], r["tipo_cluster"])] = int(r["N"])
print("\nRazon qgrad/2e10 vs T (N por T; compartido T>=1500 no se muestra):")
print(f"{'T':>5s}  {'Cu ico BOO':>18s}  {'Cu ico Voronoi':>18s}  {'Zr Z16':>18s}")
for T in (900, 800, 700, 600, 500, 400, 300):
    s = []
    for sp, tc in (("Cu", "ico_boo"), ("Cu", "ico_voronoi"), ("Zr", "Z16_voronoi")):
        a, b = N[("qgrad", T, sp, tc)], N[("2e10", T, sp, tc)]; s.append(f"{a:5d}/{b:5d}={a / b:5.2f}")
    print(f"{T:5d}  " + "  ".join(f"{x:>18s}" for x in s))
