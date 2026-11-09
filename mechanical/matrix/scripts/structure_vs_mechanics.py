#!/usr/bin/env python3
"""Vincula la estructura local (Voronoi por cuadro, de structure_evolution.py) con la respuesta mecanica (curva sigma-eps, resumen.csv, eta por atomo).
Muestras: 2e10, 1e11, 1e12, qgrad_ext; 300 K; tasas de deformacion 1e8, 1e9, 1e10; x, y, z. (qgrad excluida por pedido.)  Solo lectura; escribe CSV y figuras nuevas.
Salidas en matrix/: estructura_vs_mecanica.csv (una fila por caso) y figs/estructura_*.png
  estructura_CN_hist.png         histogramas P(CN) por especie a eps = 0 y 0.01 (estilo Jaroszewicz), 2e10 y qgrad_ext, 1e9
  estructura_dCN.png             cambio P(CN; eps) - P(CN; 0) a eps = 0.01, pico y 0.20, por muestra y especie, 1e9 (barras de error = desvio x,y,z)
  estructura_evolucion.png       sigma(eps), fraccion de ico Cu-centrado, f/f0 y <CN> global vs eps (columnas = tasa de deformacion)
  estructura_inicial_vs_mecanica.png   estado inicial (Ep, f_ico Cu, <CN> Zr, CN=15 Zr) vs sigma_max, delta_sigma, eps(sigma_max)
  estructura_perdida_vs_ablandamiento.png   ico perdido hasta el pico / hasta 0.20 vs delta_sigma; eps de maxima perdida vs eps(sigma_max)
  estructura_supervivencia_eta.png   supervivencia de los ico iniciales y eta medio de ico iniciales vs resto, por ventana de deformacion
Uso (venv_structural): python structure_vs_mechanics.py"""
import os, csv, glob, json
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); F = f"{M}/figs"
SAMPLES = ["2e10", "1e11", "1e12", "qgrad_ext"]; ER = ["1e8", "1e9", "1e10"]; DIRS = "xyz"
COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad_ext": "#b5338a"}; LS = {"qgrad_ext": "--"}
MK = {"1e8": "s", "1e9": "o", "1e10": "^"}; SPN = {1: "Zr", 2: "Cu", 3: "Ag"}; GRID = "#e6e5e1"
res = {r["case_id"]: r for r in csv.DictReader(open(f"{M}/resumen.csv"))}
def cid(s, er, d): return f"{s}_300K_{er}_{d}"
def clean(ax):
    ax.grid(color=GRID, lw=.7); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

# ---- carga ----
D = {}   # case -> dict
for s in SAMPLES:
    for er in ER:
        for d in DIRS:
            c = cid(s, er, d); p = os.path.join(M, c); f = f"{p}/analysis/structure/voro_frames.npz"
            if not os.path.exists(f): continue
            z = np.load(f); cv = np.loadtxt(glob.glob(f"{p}/*.curve.dat")[0], comments="#"); col = {"x": 2, "y": 3, "z": 4}[d]
            eta = np.load(f"{p}/analysis/localization/eta_incremental_ventanas0.02.npz")
            assert (eta["id"] == z["id"]).all()
            D[c] = dict(s=s, er=er, d=d, t=z["type"], cn=z["cn"], cls=z["cls"], eps=cv[:, 1], sig=cv[:, col], pe0=cv[:5, 6].mean(), eta=eta["eta"], eta_c=0.5 * (eta["eps_ini"] + eta["eps_fin"]), r=res[c])
print(len(D), "casos cargados")
def kpk(c): return int(round(float(D[c]["r"]["eps_sigma_max"]) / 0.005))
def fico(c, k, sp=2): m = D[c]["t"] == sp; return (D[c]["cls"][k][m] == 1).mean()
def cnm(c, k, sp=None): m = (D[c]["t"] == sp) if sp else np.ones(len(D[c]["t"]), bool); return D[c]["cn"][k][m].mean()
def pcn(c, k, sp=None):
    m = (D[c]["t"] == sp) if sp else np.ones(len(D[c]["t"]), bool); return np.bincount(D[c]["cn"][k][m], minlength=24)[6:24] / m.sum()
def ms(v): v = np.array(v, float); return v.mean(), (v.std(ddof=1) if len(v) > 1 else 0.0)

# ---- tabla por caso ----
rows = []
for c in D:
    x = D[c]; kp = kpk(c); ico0 = (x["cls"][0] == 1) & (x["t"] == 2); cu = x["t"] == 2
    f = np.array([fico(c, k) for k in range(41)]); dfdk = np.diff(f)
    ksm = np.convolve(dfdk, np.ones(3) / 3, mode="same"); kmax = int(np.argmin(ksm[:30])) + 0.5   # cuadro de maxima perdida (suavizado 3 cuadros)
    rows.append(dict(case_id=c, sample=x["s"], erate=x["er"], dir=x["d"], PE0=x["pe0"], sigma_max=float(x["r"]["sigma_max_GPa"]), delta_sigma=float(x["r"]["delta_sigma_GPa"]), E=float(x["r"]["E_GPa"]),
                     eps_peak=float(x["r"]["eps_sigma_max"]), fico_Cu_0=f[0], fico_Cu_peak=f[kp], fico_Cu_020=f[40], fico_Cu_base=f[1], loss_peak_rel=(f[0] - f[kp]) / f[0], loss_020_rel=(f[0] - f[40]) / f[0],
                     eps_maxloss=0.005 * kmax, cnZr_0=cnm(c, 0, 1), cnCu_0=cnm(c, 0, 2), cnAll_0=cnm(c, 0), cn15Zr_0=pcn(c, 0, 1)[9], cn12Cu_0=pcn(c, 0, 2)[6],
                     dcnAll_001=cnm(c, 2) - cnm(c, 0), dcnAll_peak=cnm(c, kp) - cnm(c, 0), dcnAll_020=cnm(c, 40) - cnm(c, 0),
                     surv_base=(x["cls"][1][ico0] == 1).mean(), surv_peak=(x["cls"][kp][ico0] == 1).mean(), surv_020=(x["cls"][40][ico0] == 1).mean(),
                     recr_base=(x["cls"][1][cu & ~ico0] == 1).mean(), recr_peak=(x["cls"][kp][cu & ~ico0] == 1).mean(), recr_020=(x["cls"][40][cu & ~ico0] == 1).mean()))
with open(f"{M}/estructura_vs_mecanica.csv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows([{k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()} for r in rows])
R = defaultdict(list)
for r in rows: R[(r["sample"], r["erate"])].append(r)
agg = lambda s, er, k: ms([r[k] for r in R[(s, er)]])

# ---- Fig 1: histogramas P(CN) estilo Jaroszewicz ----
CN = np.arange(6, 24); GR = [("all", None), ("Zr", 1), ("Cu", 2), ("Ag", 3)]; er = "1e9"
fg, axs = plt.subplots(2, 4, figsize=(17, 6.6), constrained_layout=True, sharex=True)
for i, s in enumerate(["2e10", "qgrad_ext"]):
    for j, (nm, sp) in enumerate(GR):
        ax = axs[i, j]
        for k, off, col, lab in ((0, -0.2, "#8a8985", "ε = 0"), (2, 0.2, "#e0a82e", "ε = 0.01")):
            h = np.mean([pcn(cid(s, er, d), k, sp) for d in DIRS], axis=0); ax.bar(CN + off, h, width=0.4, color=col, label=lab)
        ax.set_title(f"{nm} — {s}", fontsize=10.5, color=COL[s]); clean(ax)
        if j == 0: ax.set_ylabel("fracción de átomos")
        if i == 1: ax.set_xlabel("número de coordinación (Voronoi)")
        if i == 0 and j == 0: ax.legend(frameon=False, fontsize=9)
fg.suptitle("Distribución de CN a ε = 0 y 1 % (300 K, 1e9 s⁻¹; media de x, y, z). A esta escala los dos estados prácticamente coinciden: ver el cambio en estructura_dCN.png", fontsize=11)
fg.savefig(f"{F}/estructura_CN_hist.png", dpi=140); plt.close(fg)

# ---- Fig 2: dP(CN) ----
fg, axs = plt.subplots(4, 4, figsize=(17, 12), constrained_layout=True, sharex=True)
for i, s in enumerate(SAMPLES):
    kp = int(round(np.mean([kpk(cid(s, er, d)) for d in DIRS])))
    for j, (nm, sp) in enumerate(GR):
        ax = axs[i, j]; ax.axhline(0, color="#8a8985", lw=.8)
        for k, col, lab, off in ((2, "#8a8985", "ε = 0.01", -0.27), (kp, "#2aa9c9", f"ε ≈ pico ({0.005*kp:.3f})", 0.0), (40, "#0b0b0b", "ε = 0.20", 0.27)):
            dd = np.array([pcn(cid(s, er, d), k, sp) - pcn(cid(s, er, d), 0, sp) for d in DIRS]); ax.bar(CN + off, dd.mean(0), width=0.26, color=col, yerr=dd.std(0, ddof=1), error_kw=dict(lw=.7, capsize=1.5), label=lab)
        ax.set_title(f"{nm} — {s}", fontsize=10.5, color=COL[s]); clean(ax)
        if j == 0: ax.set_ylabel("P(CN; ε) − P(CN; 0)")
        if i == 3: ax.set_xlabel("número de coordinación")
        if i == 0 and j == 0: ax.legend(frameon=False, fontsize=8.5)
fg.suptitle("Cambio de la distribución de CN respecto de ε = 0 (300 K, 1e9 s⁻¹). Barras de error: desvío entre x, y, z del mismo vidrio", fontsize=11)
fg.savefig(f"{F}/estructura_dCN.png", dpi=140); plt.close(fg)

# ---- Fig 3: evolucion con eps ----
fg, axs = plt.subplots(4, 3, figsize=(16, 12), constrained_layout=True, sharex="col")
for j, er in enumerate(ER):
    for s in SAMPLES:
        cs = [cid(s, er, d) for d in DIRS if cid(s, er, d) in D]; e = D[cs[0]]["eps"]; n = min(len(D[c]["sig"]) for c in cs); sg = np.mean([D[c]["sig"][:n] for c in cs], axis=0)
        f = np.array([[fico(c, k) for k in range(41)] for c in cs]); cm = np.array([[cnm(c, k) for k in range(41)] for c in cs]); ek = 0.005 * np.arange(41)
        axs[0, j].plot(e[:n], sg, color=COL[s], ls=LS.get(s, "-"), lw=1.5, label=s)
        axs[1, j].errorbar(ek, f.mean(0), yerr=f.std(0, ddof=1), color=COL[s], ls=LS.get(s, "-"), lw=1.5, marker="o", ms=3, capsize=1.5)
        axs[2, j].plot(ek, f.mean(0) / f.mean(0)[0], color=COL[s], ls=LS.get(s, "-"), lw=1.5, marker="o", ms=3)
        axs[3, j].plot(ek, cm.mean(0) - cm.mean(0)[0], color=COL[s], ls=LS.get(s, "-"), lw=1.5, marker="o", ms=3)
        for a in axs[:, j]: a.axvline(np.mean([float(D[c]["r"]["eps_sigma_max"]) for c in cs]), color=COL[s], lw=.8, alpha=.5)
    axs[0, j].set_title(f"tasa de deformación {er} s⁻¹", fontsize=11)
    for i, yl in enumerate(["σ del eje de tracción (GPa)", "fracción de Cu con ico <0,0,12,0>", "f_ico(Cu) / f_ico(Cu; ε = 0)", "⟨CN⟩ global − ⟨CN⟩(ε = 0)"]):
        clean(axs[i, j]); axs[i, 0].set_ylabel(yl)
    axs[3, j].set_xlabel("deformación ε")
axs[0, 0].legend(frameon=False, fontsize=9)
fg.suptitle("Estructura local vs curva σ–ε (300 K). Líneas verticales finas: ε(σmax) de cada muestra. Barras: desvío entre x, y, z", fontsize=11)
fg.savefig(f"{F}/estructura_evolucion.png", dpi=140); plt.close(fg)

# ---- Fig 4: estado inicial vs mecanica ----
XS = [("PE0", "Ep inicial (eV/átomo)"), ("fico_Cu_0", "f_ico(Cu) inicial"), ("cnZr_0", "⟨CN⟩ de Zr inicial"), ("cn15Zr_0", "fracción de Zr con CN = 15 inicial")]
YS = [("sigma_max", "σmax (GPa)"), ("delta_sigma", "Δσ = σmax − σflujo (GPa)"), ("eps_peak", "ε(σmax)")]
fg, axs = plt.subplots(3, 4, figsize=(17, 10), constrained_layout=True)
corr = []
for i, (yk, yl) in enumerate(YS):
    for j, (xk, xl) in enumerate(XS):
        ax = axs[i, j]
        for s in SAMPLES:
            for er in ER:
                (mx, sx), (my, sy) = agg(s, er, xk), agg(s, er, yk); ax.errorbar(mx, my, xerr=sx, yerr=sy, color=COL[s], marker=MK[er], ms=7, mfc=COL[s] if er == "1e9" else "none", capsize=2, lw=0, elinewidth=1)
        for er in ER:
            xx = np.array([agg(s, er, xk)[0] for s in SAMPLES]); yy = np.array([agg(s, er, yk)[0] for s in SAMPLES]); r = np.corrcoef(xx, yy)[0, 1]; corr.append((yk, xk, er, r))
        r9 = [c[3] for c in corr[-3:]]; ax.text(.03, .03, "r (1e8/1e9/1e10) = " + " / ".join(f"{v:+.2f}" for v in r9), transform=ax.transAxes, va="bottom", fontsize=8.5); clean(ax)
        ax.ticklabel_format(axis='x', useOffset=False)
        if i == 2: ax.set_xlabel(xl)
        if j == 0: ax.set_ylabel(yl)
from matplotlib.lines import Line2D
axs[0, 0].legend(handles=[Line2D([], [], color=COL[s], marker="o", ls="", label=s) for s in SAMPLES] + [Line2D([], [], color="k", marker=MK[e], mfc="none", ls="", label=f"ε̇ = {e}") for e in ER], frameon=False, fontsize=8.5, loc="lower right")
fg.suptitle("Estado estructural inicial (ε = 0) vs respuesta mecánica a 300 K. Solo 4 vidrios (un punto por muestra y tasa): r es indicativo, n = 4", fontsize=11)
fg.savefig(f"{F}/estructura_inicial_vs_mecanica.png", dpi=140); plt.close(fg)

# ---- Fig 5: perdida de ico vs ablandamiento ----
fg, axs = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
for s in SAMPLES:
    for er in ER:
        for r in R[(s, er)]:
            kw = dict(color=COL[s], marker=MK[er], ms=6, mfc=COL[s] if er == "1e9" else "none", lw=0)
            axs[0].plot(r["fico_Cu_0"] - r["fico_Cu_peak"], r["delta_sigma"], **kw); axs[1].plot(r["fico_Cu_0"] - r["fico_Cu_020"], r["delta_sigma"], **kw)
            axs[2].plot(r["fico_Cu_peak"] - r["fico_Cu_020"], r["sigma_max"] - float(res[r["case_id"]]["sigma_flujo_GPa"]), **kw)
axs[0].set(xlabel="ico(Cu) perdido entre ε = 0 y el pico (fracción de los Cu)", ylabel="Δσ (GPa)"); axs[1].set(xlabel="ico(Cu) perdido entre ε = 0 y 0.20 (fracción de los Cu)", ylabel="Δσ (GPa)")
axs[2].set(xlabel="ico(Cu) perdido entre el pico y 0.20 (fracción de los Cu)", ylabel="caída post-pico σmax − σflujo (GPa)")
for a in axs: clean(a)
axs[0].legend(handles=[Line2D([], [], color=COL[s], marker="o", ls="", label=s) for s in SAMPLES] + [Line2D([], [], color="k", marker=MK[e], mfc="none", ls="", label=f"ε̇ = {e}") for e in ER], frameon=False, fontsize=8.5)
fg.suptitle("Pérdida (absoluta) de icosaedros de Cu vs ablandamiento (un punto por caso = muestra × tasa × dirección)", fontsize=11); fg.savefig(f"{F}/estructura_perdida_vs_ablandamiento.png", dpi=140); plt.close(fg)

# ---- Fig 6: supervivencia y eta ----
fg, axs = plt.subplots(2, 3, figsize=(16, 8), constrained_layout=True)
er = "1e9"
for s in SAMPLES:
    cs = [cid(s, er, d) for d in DIRS]; ek = 0.005 * np.arange(41)
    S = np.array([[((D[c]["cls"][k] == 1)[(D[c]["cls"][0] == 1) & (D[c]["t"] == 2)]).mean() for k in range(41)] for c in cs])
    Rr = np.array([[((D[c]["cls"][k] == 1)[(D[c]["cls"][0] != 1) & (D[c]["t"] == 2)]).mean() for k in range(41)] for c in cs])
    axs[0, 0].plot(ek[1:], S.mean(0)[1:], color=COL[s], ls=LS.get(s, "-"), lw=1.5, marker="o", ms=3, label=s); axs[0, 1].plot(ek[1:], Rr.mean(0)[1:], color=COL[s], ls=LS.get(s, "-"), lw=1.5, marker="o", ms=3)
    axs[0, 2].plot(ek[1:], S.mean(0)[1:] / S.mean(0)[1], color=COL[s], ls=LS.get(s, "-"), lw=1.5, marker="o", ms=3)
    for j, (nm, fn) in enumerate([("ico inicial / Cu no-ico inicial", lambda c, w: D[c]["eta"][w][((D[c]["cls"][0] == 1) & (D[c]["t"] == 2))].mean() / D[c]["eta"][w][((D[c]["cls"][0] != 1) & (D[c]["t"] == 2))].mean()),
                                  ("Cu con CN = 12 inicial / Cu resto", lambda c, w: D[c]["eta"][w][(D[c]["cn"][0] == 12) & (D[c]["t"] == 2)].mean() / D[c]["eta"][w][(D[c]["cn"][0] != 12) & (D[c]["t"] == 2)].mean()),
                                  ("Zr con CN = 15 inicial / Zr resto", lambda c, w: D[c]["eta"][w][(D[c]["cn"][0] == 15) & (D[c]["t"] == 1)].mean() / D[c]["eta"][w][(D[c]["cn"][0] != 15) & (D[c]["t"] == 1)].mean())]):
        v = np.array([[fn(c, w) for w in range(10)] for c in cs]); ec = D[cs[0]]["eta_c"]
        axs[1, j].errorbar(ec, v.mean(0), yerr=v.std(0, ddof=1), color=COL[s], ls=LS.get(s, "-"), lw=1.5, marker="o", ms=4, capsize=2); axs[1, j].axhline(1, color="#8a8985", lw=.8); axs[1, j].set_title(f"η medio: {nm}", fontsize=10); axs[1, j].set_xlabel("ε central de la ventana (Δε = 0.02)")
axs[0, 0].set(ylabel="ico(Cu) iniciales que siguen siéndolo", xlabel="ε"); axs[0, 1].set(ylabel="Cu no-ico iniciales que pasan a ico", xlabel="ε"); axs[0, 2].set(ylabel="supervivencia normalizada a ε = 0.005", xlabel="ε")
axs[1, 0].set_ylabel("cociente de η medio (ventana)")
for a in axs.flat: clean(a)
axs[0, 0].legend(frameon=False, fontsize=9)
fg.suptitle("Destino de los icosaedros iniciales y deformación local posterior (300 K, 1e9 s⁻¹; media de x, y, z). η = deformación de von Mises incremental por ventana", fontsize=11)
fg.savefig(f"{F}/estructura_supervivencia_eta.png", dpi=140); plt.close(fg)

# ---- resumen numerico ----
print("== f_ico(Cu) (media sobre x,y,z) a eps=0, en el pico, a 0.20; perdida relativa al pico y a 0.20; supervivencia a 0.20; <CN> global cambio 0.01/pico/0.20")
for s in SAMPLES:
    for er in ER:
        print(f"  {s:9s} {er:5s} f0={agg(s,er,'fico_Cu_0')[0]:.4f} f_base={agg(s,er,'fico_Cu_base')[0]:.4f} f_pk={agg(s,er,'fico_Cu_peak')[0]:.4f} f_020={agg(s,er,'fico_Cu_020')[0]:.4f} | loss_pk={agg(s,er,'loss_peak_rel')[0]*100:.1f}% loss_020={agg(s,er,'loss_020_rel')[0]*100:.1f}% | "
              f"surv_base={agg(s,er,'surv_base')[0]:.3f} surv_020={agg(s,er,'surv_020')[0]:.3f} | dCN 0.01={agg(s,er,'dcnAll_001')[0]:+.4f} pk={agg(s,er,'dcnAll_peak')[0]:+.4f} 020={agg(s,er,'dcnAll_020')[0]:+.4f} | sigma_max={agg(s,er,'sigma_max')[0]:.3f} d_sigma={agg(s,er,'delta_sigma')[0]:.3f}")
print("== correlaciones (n = 4 muestras)")
for c in corr: print(f"  {c[0]:12s} vs {c[1]:10s} @ {c[2]:5s}: r = {c[3]:+.2f}")
