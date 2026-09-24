#!/usr/bin/env python3
"""Etapa 4 / Tarea 2 -- BADF P(alpha) (Fig C).

Configs : ../../<tasa>/snaps/17.dump (glass 300 K) y 08.dump (liquido 1200 K); posiciones
          wrapped de la Etapa 0. T = 1900 - 100*(NN-1).
Cutoff  : primer minimo de g_total(r) de ESA config (../rdf/rdf_all.csv), aplicado a todas
          las especies (Zr/Cu/Ag) como vecinos.
Angulos : para cada atomo central i, todos los pares (j,k) de vecinos dentro de rc, con
          vectores r_j - r_i en minima imagen (caja cubica periodica). Variantes: center=all
          (todo i) y center=Cu (solo type 2 = Cu; orden EAM Zr Cu Ag).
Normalizacion: histograma de 1 deg (180 bins en [0,180]) de TODOS los angulos acumulados;
          P(alpha_k) = 180 * cuentas_k / N_total  (alpha en grados = centro de bin)
          => densidad respecto de x = alpha/pi: integral_0^1 P dx = sum_k P_k/180 = 1, es decir
          MEDIA de P sobre [0,180 deg] = 1 y distribucion uniforme en alpha -> P = 1.
          Es la escala de la Fig 4 de Ryltsev (vidrio ~3.3 a 60 deg); P_por_radian = P/pi.
          Sin peso 1/sin(alpha). Misma funcion para glass, liquido e icosaedro/hcp ideales.
Salida  : badf.csv (rate,state,center,alpha,P), badf_refs.csv (ref,alpha,P), badf_meta.csv
          ../fig_C_badf.png
Sin pandas: numpy + scipy + csv.
"""
import csv
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.signal import savgol_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ETAPA4 = HERE.parent
SA = ETAPA4.parent                                  # structural_analysis/
RATES = ["1e12", "1e11", "2e10"]
STATES = {"glass300": (17, 300), "liquid": (8, 1200)}   # (NN, T)
TYPE_CU = 2
NB = 180
EDGES = np.linspace(0.0, 180.0, NB + 1)
ALPHA = 0.5 * (EDGES[:-1] + EDGES[1:])
DX = 1.0 / NB                                       # ancho de bin en x = alpha/pi
RATE_LABEL = {"1e12": r"$10^{12}$ K/s", "1e11": r"$10^{11}$ K/s", "2e10": r"$2\times10^{10}$ K/s"}
RATE_COLOR = {"1e12": "#d62728", "1e11": "#e6b800", "2e10": "#1f3fd6"}


# ------------------------------------------------------------------ io
def read_dump(path):
    with open(path) as f:
        for _ in range(3):
            f.readline()
        n = int(f.readline())
        f.readline()
        lo, hi = [], []
        for _ in range(3):
            a, b = f.readline().split()[:2]
            lo.append(float(a)); hi.append(float(b))
        f.readline()
        arr = np.loadtxt(f)
    assert arr.shape == (n, 5)
    arr = arr[np.argsort(arr[:, 0])]
    L = np.array(hi) - np.array(lo)
    assert np.allclose(L, L[0], atol=1e-9), "caja no cubica"
    return arr[:, 1].astype(int), arr[:, 2:5] - np.array(lo), float(L[0])


def load_rdf():
    d = {}
    with open(ETAPA4 / "rdf" / "rdf_all.csv") as f:
        for row in csv.DictReader(f):
            d.setdefault((row["rate"], int(row["T"])), []).append((float(row["r"]), float(row["g_total"])))
    return {k: np.array(v) for k, v in d.items()}


def load_qc_L():
    with open(ETAPA4 / "rdf" / "rdf_qc.csv") as f:
        return {(r["rate"], int(r["T"])): float(r["L"]) for r in csv.DictReader(f)}


# ------------------------------------------------------------------ primer minimo
def first_min(r, g):
    """Primer minimo local de g(r) tras el primer pico. Suavizado Savitzky-Golay (9 bins = 0.18 A,
    polinomio 3); pico = maximo del suavizado en [2.0, 3.6] A; minimo = primer punto i > pico
    (r > r_pico + 0.2 A) con gs[i-1] > gs[i] <= gs[i+1]. Devuelve (r_pico, r_min, g_min_suavizado)."""
    gs = savgol_filter(g, 9, 3)
    win = (r >= 2.0) & (r <= 3.6)
    ipk = np.where(win)[0][np.argmax(gs[win])]
    for i in range(ipk + 1, len(r) - 1):
        if r[i] > r[ipk] + 0.2 and gs[i - 1] > gs[i] <= gs[i + 1]:
            return r[ipk], r[i], gs[i]
    raise RuntimeError("no se encontro primer minimo")


# ------------------------------------------------------------------ BADF
def angles_from_vectors(V):
    """V: (n,3) vectores a los n vecinos. Devuelve los n(n-1)/2 angulos (grados) entre pares."""
    u = V / np.linalg.norm(V, axis=1, keepdims=True)
    c = np.clip(u @ u.T, -1.0, 1.0)
    iu = np.triu_indices(len(V), k=1)
    return np.degrees(np.arccos(c[iu]))


def normalize(counts):
    return counts / (counts.sum() * DX)


def badf_config(types, pos, L, rc):
    tree = cKDTree(pos, boxsize=L)
    nbrs = tree.query_ball_point(pos, r=rc)
    cnt = {"all": np.zeros(NB), "Cu": np.zeros(NB)}
    ncen = {"all": 0, "Cu": 0}
    cn_sum = {"all": 0, "Cu": 0}
    for i, nb in enumerate(nbrs):
        nb = [j for j in nb if j != i]
        if len(nb) < 2:
            continue
        d = pos[nb] - pos[i]
        d -= L * np.round(d / L)                      # minima imagen
        h, _ = np.histogram(angles_from_vectors(d), bins=EDGES)
        cnt["all"] += h; ncen["all"] += 1; cn_sum["all"] += len(nb)
        if types[i] == TYPE_CU:
            cnt["Cu"] += h; ncen["Cu"] += 1; cn_sum["Cu"] += len(nb)
    meta = {c: dict(n_centers=ncen[c], mean_CN=cn_sum[c] / ncen[c], n_angles=int(cnt[c].sum())) for c in cnt}
    return {c: normalize(cnt[c]) for c in cnt}, meta


# ------------------------------------------------------------------ referencias ideales
def ico_neighbors():
    phi = (1 + 5 ** 0.5) / 2
    v = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            v += [(0, s1, s2 * phi), (s1, s2 * phi, 0), (s2 * phi, 0, s1)]
    v = np.array(v, float)
    assert len(v) == 12
    return v


def hcp_neighbors():
    a = 1.0
    v = [(a * np.cos(np.radians(60 * k)), a * np.sin(np.radians(60 * k)), 0.0) for k in range(6)]
    zc = a * np.sqrt(2.0 / 3.0)                       # c/2 con c/a = sqrt(8/3)
    rho = a / np.sqrt(3.0)
    v += [(rho * np.cos(np.radians(t)), rho * np.sin(np.radians(t)), +zc) for t in (30, 150, 270)]
    v += [(rho * np.cos(np.radians(t)), rho * np.sin(np.radians(t)), -zc) for t in (90, 210, 330)]
    v = np.array(v)
    assert np.allclose(np.linalg.norm(v, axis=1), a)   # 12 vecinos a distancia a
    return v


def ref_badf(V):
    h, _ = np.histogram(angles_from_vectors(V), bins=EDGES)
    return normalize(h.astype(float))


def gauss_smooth(P, sigma_deg=3.0):
    """Solo para el PNG: convolucion gaussiana (sigma en grados, bins de 1 deg), area conservada."""
    k = np.arange(-4 * int(sigma_deg) - 1, 4 * int(sigma_deg) + 2)
    w = np.exp(-0.5 * (k / sigma_deg) ** 2); w /= w.sum()
    return np.convolve(np.pad(P, len(k), mode="reflect"), w, mode="same")[len(k):-len(k)]


# ------------------------------------------------------------------ main
rdf = load_rdf()
qcL = load_qc_L()
rows, meta_rows, res = [], [], {}
for rate in RATES:
    for state, (nn, T) in STATES.items():
        types, pos, L = read_dump(SA / rate / "snaps" / f"{nn:02d}.dump")
        assert abs(L - qcL[(rate, T)]) < 1e-3, f"caja dump {L} vs rdf_qc {qcL[(rate, T)]}"
        a = rdf[(rate, T)]
        r_pk, rc, g_rc = first_min(a[:, 0], a[:, 1])
        P, meta = badf_config(types, pos, L, rc)
        for c in ("all", "Cu"):
            res[(rate, state, c)] = P[c]
            for al, p in zip(ALPHA, P[c]):
                rows.append([rate, state, c, f"{al:.1f}", f"{p:.6f}"])
            meta_rows.append([rate, state, c, T, f"{L:.4f}", f"{r_pk:.3f}", f"{rc:.3f}", f"{g_rc:.4f}",
                              meta[c]["n_centers"], f"{meta[c]['mean_CN']:.3f}", meta[c]["n_angles"],
                              f"{P[c].sum() * DX:.6f}"])
            print(f"{rate:5s} {state:9s} {c:3s} T={T:4d} r_pico={r_pk:.2f} rc={rc:.2f} g(rc)={g_rc:.3f} "
                  f"centros={meta[c]['n_centers']:5d} <CN>={meta[c]['mean_CN']:.2f} area={P[c].sum()*DX:.4f}")

with open(HERE / "badf.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rate", "state", "center", "alpha", "P"]); w.writerows(rows)
with open(HERE / "badf_meta.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rate", "state", "center", "T", "L", "r_pico", "rc", "g_rc", "n_centers", "mean_CN",
                "n_angles", "area"])
    w.writerows(meta_rows)

refs = {"ico": ref_badf(ico_neighbors()), "hcp": ref_badf(hcp_neighbors())}
with open(HERE / "badf_refs.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["ref", "alpha", "P"])
    for k, P in refs.items():
        for al, p in zip(ALPHA, P):
            w.writerow([k, f"{al:.1f}", f"{p:.6f}"])
for k, P in refs.items():
    top = ALPHA[np.argsort(P)[::-1][:4]]
    print(f"ref {k}: area={P.sum()*DX:.4f}  bins mas altos (deg)={sorted(top)}")

# ------------------------------------------------------------------ Fig C
fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.6), sharey=True)
titles = {"all": "all atoms as centre", "Cu": "Cu-centred"}
for ax, c in zip(axes, ("all", "Cu")):
    for i, rate in enumerate(RATES):
        ax.plot(ALPHA, res[(rate, "liquid", c)], color="0.15", lw=1.0, ls="-", alpha=0.55,
                label="liquid 1200 K (3 rates)" if i == 0 else None)
    for rate in RATES:
        ax.plot(ALPHA, res[(rate, "glass300", c)], color=RATE_COLOR[rate], lw=1.8,
                label=f"glass 300 K, {RATE_LABEL[rate]}")
    ax.plot(ALPHA, gauss_smooth(refs["ico"]), color="olive", ls="--", lw=1.4, label="ico (ideal, σ=3° blur)")
    ax.plot(ALPHA, gauss_smooth(refs["hcp"]), color="magenta", ls="--", lw=1.4, label="hcp (ideal, σ=3° blur)")
    ax.set_xlim(0, 180); ax.set_ylim(0, 4)
    ax.set_xticks(range(0, 181, 60))
    ax.set_xlabel(r"$\alpha$, deg")
    ax.set_title(titles[c], fontsize=10)
axes[0].set_ylabel(r"$P(\alpha)$")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0.12, 1, 1))
fig.savefig(ETAPA4 / "fig_C_badf.png", dpi=200)
print("fig_C escrita")
