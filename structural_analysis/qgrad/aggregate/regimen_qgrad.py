#!/usr/bin/env python3
"""Etapa 3 (paso 2): razon por regimen qgrad / referencia (2e10 = tronco, mas 1e11 y 1e12 de contexto).
MISMO rigor que step9 (no se cambia ningun criterio):
  - filtro N>=20 POR T en AMBAS curvas, T que falla se descarta; si ninguno pasa -> no calculable
  - razon = N_pool(qgrad)/N_pool(ref) (>1: mas orden en qgrad); error = razon*sqrt(1/N1+1/N2) (piso Poisson)
  - ventana compartida = interseccion exacta: liquido = {T > Tg_a+100} & {T > Tg_b+100}; vitreo = {T < Tg_a-100} & {T < Tg_b-100}
  - Tg = media de Tg(V) y Tg(Ep) (misma convencion que step6/step9). Tg de qgrad = Tg del PERFIL (sin tasa unica).
Adicion propia de qgrad: los T>=1500 K de qgrad son la trayectoria 2e10 (mismas configs, no independientes) -> se EXCLUYEN
de todos los pools (aparecen en 'excluidos_compartidos'). En vitreo, 300 K de qgrad = 17r si existe, si no 17u (sin relajar) y
se anota en la fila. Uso: python regimen_qgrad.py [--300K 17u|17r]"""
import sys, csv, math
NMIN, BUF = 20, 100.0
opt300 = sys.argv[sys.argv.index("--300K") + 1] if "--300K" in sys.argv else None
Tg = {"2e10": (646.1 + 637.2) / 2, "1e11": (669.3 + 663.0) / 2, "1e12": (703.2 + 697.9) / 2}
std = [r for r in csv.DictReader(open("../thermo/Tg_qgrad_y_sensibilidad.csv")) if r["curva"] == "qgrad" and r["ventanas"].startswith("ventanas estandar")][0]
Tg["qgrad"] = (float(std["Tg_V_K"]) + float(std["Tg_Ep_K"])) / 2
print("Tg usadas:", {k: round(v, 2) for k, v in Tg.items()})

N, comp, cfg300 = {}, {}, set()
for r in csv.DictReader(open("tabla_compare_rates_4curvas.csv")):
    T = int(r["T"]); key = (r["rate"], T, r["especie"], r["tipo_cluster"])
    if r["rate"] == "qgrad" and T == 300:            # variantes 17u / 17r
        cfg300.add(r["config_qgrad"]); key = (r["rate"] + ":" + r["config_qgrad"], T, r["especie"], r["tipo_cluster"])
    N[key] = int(r["N"]); comp[key] = int(r["compartido_con_2e10"])
lab300 = opt300 or ("17r" if "17r" in cfg300 else "17u")
print("300 K de qgrad =", lab300, "(disponibles:", sorted(cfg300), ")")
Ts = sorted({k[1] for k in N}, reverse=True)
combos = sorted({(k[2], k[3]) for k in N})

def n_of(rate, T, sp, tc):
    if rate == "qgrad" and T == 300: rate = "qgrad:" + lab300
    return N.get((rate, T, sp, tc))

out = []
for ref in ["2e10", "1e11", "1e12"]:
    win = {"liquido": [T for T in Ts if T > Tg["qgrad"] + BUF and T > Tg[ref] + BUF],
           "vitreo": [T for T in Ts if T < Tg["qgrad"] - BUF and T < Tg[ref] - BUF]}
    for reg, wT in win.items():
        excl = [T for T in wT if T >= 1500]           # T compartidos con 2e10
        use_T = [T for T in wT if T < 1500]
        for sp, tc in combos:
            per = [(T, n_of("qgrad", T, sp, tc), n_of(ref, T, sp, tc)) for T in use_T]
            per = [(T, a, b) for T, a, b in per if a is not None and b is not None]
            kept = [(T, a, b) for T, a, b in per if a >= NMIN and b >= NMIN]
            drop = [T for T, a, b in per if not (a >= NMIN and b >= NMIN)]
            n1, n2 = sum(a for _, a, _ in kept), sum(b for _, _, b in kept)
            if n1 and n2:
                r_ = n1 / n2; e_ = r_ * math.sqrt(1 / n1 + 1 / n2); s_ = (r_ - 1) / e_
                vals = (n1, n2, f"{r_:.4f}", f"{e_:.4f}", f"{s_:.2f}")
            else:
                vals = ("", "", "", "", "")
            nota = ("" if kept else "NO CALCULABLE (N<%d en todos los T)" % NMIN)
            if reg == "vitreo" and 300 in use_T: nota = (nota + "; " if nota else "") + f"300K qgrad={lab300}"
            out.append((ref, sp, tc, reg, f"{min(use_T)}-{max(use_T)}" if use_T else "", " ".join(str(t) for t, _, _ in kept),
                        " ".join(map(str, drop)), " ".join(map(str, excl)), *vals, nota))
hdr = ["referencia", "especie", "tipo_cluster", "regimen", "ventana_T_nuevos", "T_usados", "T_descartados_N<20",
       "excluidos_compartidos", "N_pool_qgrad", "N_pool_ref", "razon_qgrad/ref", "error", "(razon-1)/error", "nota"]
with open("tabla_regimen_qgrad.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(out)
for ref in ["2e10"]:
    print(f"\n=== qgrad / {ref} ===")
    for o in out:
        if o[0] == ref:
            print(f"{o[1]:2s} {o[2]:12s} {o[3]:8s} T[{o[5]}] pool {o[8]}/{o[9]}  razon {o[10]} +- {o[11]}  ({o[12]}σ) {o[13]}")
