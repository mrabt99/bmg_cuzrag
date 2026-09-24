#!/usr/bin/env python3
"""Tarea Q -- aplica el filtro N>=20 POR T (ambas tasas) a la fila Cu-BOO liquido
(ventana compartida 1900-900 K), que venia sin filtrar del step6 original.
Solo modifica esa fila de tabla_regimen_pool.csv; el resto queda igual."""
import csv, math
NMIN = 20
rows = list(csv.DictReader(open("tabla_clusters_51.csv")))
Ts = list(range(1900, 800, -100))
cnt = {t: {r["tasa"]: int(r["N_ico_BOO"]) for r in rows
           if r["especie"] == "Cu" and int(r["T_nominal"]) == t and r["tasa"] in ("2e10", "1e12")} for t in Ts}
kept = [t for t in Ts if cnt[t]["2e10"] >= NMIN and cnt[t]["1e12"] >= NMIN]
n1, n2 = sum(cnt[t]["2e10"] for t in kept), sum(cnt[t]["1e12"] for t in kept)
print("desglose 2e10/1e12:", {t: f"{cnt[t]['2e10']}/{cnt[t]['1e12']}" for t in Ts})
print(f"pool sin filtro {sum(cnt[t]['2e10'] for t in Ts)}/{sum(cnt[t]['1e12'] for t in Ts)}; sobreviven {len(kept)}/{len(Ts)}: {kept}")
assert kept, "ningun T pasa -> habria que marcar no calculable"
r = n1 / n2; e = r * math.sqrt(1 / n1 + 1 / n2)
print(f"razon {r:.4f} +- {e:.4f}  (N {n1}/{n2})")

nota = (f"BOO: filtro N>={NMIN} por T deja 1 de 11 T ({' '.join(map(str, kept))} K; el resto <{NMIN} en alguna tasa: "
        "0-3 atomos a 1900-1200 K, 2-15 a 1100-1000 K). Razon de UN solo T, en el borde frio de la ventana. "
        "Lectura fisica: el orden icosaedrico es un fenomeno de liquido sobreenfriado que se construye al acercarse a Tg; "
        "en la ventana liquida, mayormente lejos de Tg, se espera poca poblacion, no es solo un problema de conteo "
        "(la poblacion BOO crece al bajar T: 0 a 1900-1700 K -> 36 a 900 K). El criterio BOO es ademas mas estricto "
        "que Voronoi <0,0,12,0>, que ya da 107 atomos a 1900 K. Valor previo sin filtro: 1.2182+-0.2217 (67/55, ~5 por T)")
tab = list(csv.reader(open("tabla_regimen_pool.csv", newline="")))
h = tab[0]; idx = {c: i for i, c in enumerate(h)}
n_mod = 0
for row in tab[1:]:
    if row[idx["especie"]] == "Cu" and row[idx["metodo"]] == "BOO" and row[idx["regimen"]] == "liquido":
        row[idx["N_pool_2e10"]], row[idx["N_pool_1e12"]] = str(n1), str(n2)
        row[idx["razon"]], row[idx["error"]] = f"{r:.4f}", f"{e:.4f}"
        row[idx["nota"]] = nota; n_mod += 1
assert n_mod == 1
csv.writer(open("tabla_regimen_pool.csv", "w", newline="")).writerows(tab)
print("-> tabla_regimen_pool.csv: 1 fila modificada (Cu,BOO,liquido)")
