"""Cierre de la Fase 1: chequeos pendientes sobre las series de cal_1200K_<tasa>.
Uso: venv_structural/bin/python fase1_cierre.py [tasa ...]     (por defecto 2e10 1e11 1e12)
Revisa: (a) media de pxy/pxz/pyz y de J (impacto de no restar la media en acf),
        (b) compatibilidad de las tres componentes de cizalla en eta,
        (c) tiempo de decaimiento de la ACF y ventana, (d) replicas necesarias.
"""
import sys, numpy as np

T = 1200.0
kB = 1.380649e-23; eV = 1.602176634e-19
NBLK = 10
TAGS = ["cal_r1", "cal_r2", "cal_r3"]
rates = sys.argv[1:] or ["2e10", "1e11", "1e12"]


def acf(x, nlag, sub_mean=False):
    n = len(x)
    if sub_mean:
        x = x - x.mean()
    f = np.fft.rfft(x, 2 * n)
    return np.fft.irfft(f * np.conj(f))[:nlag] / (n - np.arange(nlag))


def blocks(cols, nlag, sub_mean=False):
    n = len(cols[0]); m = n // NBLK
    return np.array([np.mean([acf(c[b*m:(b+1)*m], nlag, sub_mean) for c in cols], axis=0)
                     for b in range(NBLK)])


for rate in rates:
    d = f"cal_1200K_{rate}"
    V = float(open(f"{d}/cal_vmean.inc").read().split()[3])
    print(f"\n################ {rate}   V={V:.6g} A^3")

    P, J = {}, {}
    for tg in TAGS:
        P[tg] = np.loadtxt(f"{d}/{tg}_press.dat", comments="#")
        J[tg] = np.loadtxt(f"{d}/{tg}_hflux.dat", comments="#")
    dtp = (P[TAGS[0]][1, 0] - P[TAGS[0]][0, 0]) * 1e-3
    dtj = (J[TAGS[0]][1, 0] - J[TAGS[0]][0, 0]) * 1e-3
    nlp, nlj = int(20 / dtp), int(20 / dtj)
    tp, tj = np.arange(nlp) * dtp, np.arange(nlj) * dtj
    pref = V * 1e-30 / (kB * T) * 1e-12 * dtp * 1e3           # -> mPa s
    kpref = 1.0 / (V * 1e-30 * kB * T**2) * 1e-12 * dtj       # J.J = suma de 3 comps

    # ---- (a) medias ----
    print("\n(a) media / desvio de las series (media/desvio chico => restar la media no importa)")
    names = ["pxy", "pxz", "pyz"]
    for i, nm in enumerate(names):
        r = [P[tg][:, 4 + i].mean() / P[tg][:, 4 + i].std() for tg in TAGS]
        print(f"  {nm}: media/std por replica = " + "  ".join(f"{x:+.4f}" for x in r)
              + f"   (std tipica {P[TAGS[0]][:, 4+i].std():.0f} bar)")
    for i, nm in enumerate(["Jx", "Jy", "Jz"]):
        r = [J[tg][:, 1 + i].mean() / J[tg][:, 1 + i].std() for tg in TAGS]
        print(f"  {nm}: media/std por replica = " + "  ".join(f"{x:+.4f}" for x in r))

    # ---- (b) eta por componente, con y sin restar la media ----
    print("\n(b) eta (mPa s) por componente, global sobre 30 bloques  [t=3, 5, 10 ps]")
    ETA = {}
    for nm, cols in [("pxy", [4]), ("pxz", [5]), ("pyz", [6]), ("3 comp", [4, 5, 6])]:
        allb = []
        for tg in TAGS:
            B = blocks([P[tg][:, c] * 1e5 for c in cols], nlp)
            allb.append((np.cumsum(B, 1) - 0.5*B[:, :1]) * pref)
        ETA[nm] = np.concatenate(allb)                         # (30, lag)
        out = []
        for tw in (3, 5, 10):
            g = ETA[nm][:, int(tw / dtp)]
            out.append(f"{g.mean():.2f}+-{g.std(ddof=1)/np.sqrt(len(g)):.2f}")
        print(f"  {nm:7s}: " + "   ".join(out))
    allb = []
    for tg in TAGS:
        B = blocks([P[tg][:, c] * 1e5 for c in (4, 5, 6)], nlp, sub_mean=True)
        allb.append((np.cumsum(B, 1) - 0.5*B[:, :1]) * pref)
    E0 = np.concatenate(allb)
    g0, g1 = E0[:, int(10 / dtp)], ETA["3 comp"][:, int(10 / dtp)]
    print(f"  3 comp restando media (t=10 ps): {g0.mean():.3f}  vs sin restar: {g1.mean():.3f}"
          f"  -> dif {100*(g0.mean()-g1.mean())/g1.mean():+.2f} %")

    # ---- (c) decaimiento de la ACF y meseta ----
    Bm = ETA  # solo para reutilizar nombres
    Bp = np.mean([blocks([P[tg][:, c] * 1e5 for c in (4, 5, 6)], nlp).mean(0) for tg in TAGS], axis=0)
    Bp = Bp / Bp[0]
    z = np.where(Bp < 0)[0]
    e = np.where(Bp < np.exp(-1))[0]
    print(f"\n(c) ACF de presion (3 comp): cae a 1/e en t={tp[e[0]]:.3f} ps; primer cruce por 0 en "
          f"{'t=%.3f ps' % tp[z[0]] if len(z) else 'ninguno'}")
    Jb = np.mean([blocks([J[tg][:, 1 + i] * eV * 1e-10 / 1e-12 for i in range(3)], nlj).mean(0)
                  for tg in TAGS], axis=0)
    Jn = Jb / Jb[0]
    ej = np.where(Jn < np.exp(-1))[0]; zj = np.where(Jn < 0)[0]
    print(f"    ACF de J: 1/e en t={tj[ej[0]]:.3f} ps; primer cruce por 0 en "
          f"{'t=%.3f ps' % tj[zj[0]] if len(zj) else 'ninguno'}")
    print("    eta(t) global [mPa s]:  " + "  ".join(
        f"t={tw}:{ETA['3 comp'][:, int(tw/dtp)].mean():.2f}" for tw in (0.5, 1, 2, 3, 5, 7, 10, 15, 19)))
    KJ = []
    for tg in TAGS:
        B = blocks([J[tg][:, 1 + i] * eV * 1e-10 / 1e-12 for i in range(3)], nlj)
        KJ.append((np.cumsum(B, 1) - 0.5*B[:, :1]) * 3 * kpref / 3)             # acf ya promedia las 3 comps
    KJ = np.concatenate(KJ)
    print("    kappa(t) global [W/m/K]: " + "  ".join(
        f"t={tw}:{KJ[:, int(tw/dtj)].mean():.3f}" for tw in (0.5, 1, 2, 3, 5, 7, 10, 15, 19)))

    # ---- (d) replicas necesarias ----
    print("\n(d) error relativo y bloques (de 100 ps) necesarios para 3 % y 5 %")
    for nm, A, dt in [("eta", ETA["3 comp"], dtp), ("kappa", KJ, dtj)]:
        for tw in (3, 5, 10):
            g = A[:, int(tw / dt)]
            rel = g.std(ddof=1) / abs(g.mean())                # error relativo por bloque
            n3, n5 = (rel / 0.03) ** 2, (rel / 0.05) ** 2
            print(f"  {nm:5s} t={tw:2d} ps: sem={100*rel/np.sqrt(len(g)):.1f} %   "
                  f"bloques para 5 %: {n5:.0f}   para 3 %: {n3:.0f}   (hoy 30)")
