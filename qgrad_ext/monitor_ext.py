#!/usr/bin/env python3
"""Monitor de la rama qgrad_ext (SIN analisis estructural). Solo LEE logs; no lanza nada. Escribe en monitor/.
Por segmento: (1) continuidad en el borde (T, Lx, rho, Ep) entre el 'terminado' del tramo N y la 1a fila de thermo del tramo N+1 (estado exacto
del restart) y entre promedios de ventana de log.dat; (2) stderr limpio (filtra el ruido PMIx/munge/xalt) y sin halt/ERROR/lost atoms/nan en el log;
(3) check_crystallization.py (el original de qgrad) sobre cada log.dat; (4) figura Ep(T) de la rama sobre 2e10, 1e11, 1e12 y qgrad + Ep-Ep(2e10) + pendiente local,
y la ALARMA 'rama ya lineal y paralela al vidrio de 2e10' (SOLO informa; la decision de cortar es del usuario).
Uso:  python3 monitor_ext.py                                  # rama qgrad_ext: tramos 09 10 11 12 y relax que existan
      python3 monitor_ext.py --dir <carpeta> --prefix gradmin --segs 08 09 10 11 12 --start-log <log_del_tramo_previo>   # prueba sobre qgrad
Requiere numpy y matplotlib (usar venv_structural: /scratch/negishi/aaabt/venv_structural/bin/python)."""
import argparse, glob, os, re, subprocess, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); SA = "/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis"
CHECK = "/scratch/negishi/aaabt/bmg_cuzrag/rate_gradiente_gammamin/check_crystallization.py"
QGRAD_DIR = "/scratch/negishi/aaabt/bmg_cuzrag/rate_gradiente_gammamin"
NOISE = re.compile(r"PMIx|munge|psec|Host:|Framework|Component|component|installed|loaded|libraries|used on your|not found|xalt|purge|module|^-+$|^\s*$|Use \"module")
TOL = dict(T=3.0, Lx=0.013, rho=0.002, Ep=0.3e-3)      # tolerancias de la auditoria de qgrad (promedios de ventana; Ep en eV/at)
ap = argparse.ArgumentParser(); ap.add_argument("--dir", default=HERE.rsplit("/monitor", 1)[0]); ap.add_argument("--prefix", default="gradext")
ap.add_argument("--segs", nargs="*", default=["09", "10a", "10", "11", "12"]); ap.add_argument("--start-log", default=None, help="log del tramo PREVIO al primero (para su 'terminado')")
ap.add_argument("--out", default=os.path.join(HERE, "monitor")); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
D, P = a.dir, a.prefix
def logfile(seg): g = sorted(glob.glob(f"{D}/log.seg{seg}.*")); return g[-1] if g else None
def datfile(seg): p = f"{D}/quench_N65536_{P}.log.dat.seg{seg}"; return p if os.path.exists(p) else None
def read_dat(p):
    r = [l.split() for l in open(p) if l[0] != "#" and len(l.split()) == 6]; return np.array(r, float)   # step T Ep Lx rho P
def terminado(logp):
    t = open(logp).read(); m = re.findall(r"### tramo .*terminado : T=(\S+) K\s+P=(\S+) bar\s+Ep=(\S+) eV\s+Lx=(\S+) A\s+rho=(\S+) g/cc", t)
    return None if not m else dict(zip(("T", "P", "Ep", "Lx", "rho"), map(float, m[-1])))
def first_thermo(logp):
    t = open(logp).read(); m = re.search(r"Step\s+Temp\s+Press\s+Volume\s+Density\s+PotEng\s+TotEng\s*\n\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", t)
    if not m: return None
    s, T, Pp, V, rho, pe = m.groups(); return dict(step=int(s), T=float(T), P=float(Pp), Lx=float(V) ** (1 / 3), rho=float(rho), Ep=float(pe) / 65536)
print(f"== monitor_ext: dir={D} prefix={P} tramos={a.segs} ==")
rows, ok_all = [], True
prev_end = terminado(a.start_log) if a.start_log else None
prev_dat_last = None
for seg in a.segs:
    lp, dp = logfile(seg), datfile(seg)
    if not lp or not dp: print(f"tramo {seg}: (todavia no existe log/log.dat)"); break
    end = terminado(lp); ft = first_thermo(lp); dat = read_dat(dp); status = []
    print(f"\n--- tramo {seg}: {len(dat)} filas de log.dat, T {dat[0,1]:.1f} -> {dat[-1,1]:.1f} K")
    if end is None: print("   AVISO: el tramo NO termino (falta la linea '### tramo ... terminado'): en curso o fallo"); ok_all = False; status.append("sin_terminar")
    if prev_end is not None and ft is not None:                                    # borde exacto restart -> primera fila de thermo
        d = {k: ft[k] - prev_end[k] for k in ("T", "P", "Lx", "rho", "Ep")}
        good = all(abs(v) < 1e-3 * max(1, abs(prev_end[k])) or abs(v) < 1e-6 for k, v in d.items())
        print(f"   borde EXACTO (terminado tramo previo -> paso 0 de este): dT={d['T']:+.4f} dP={d['P']:+.3f} dLx={d['Lx']:+.2e} drho={d['rho']:+.2e} dEp={d['Ep']:+.2e}  {'OK' if good else 'REVISAR'}")
        ok_all &= good; status.append("borde_exacto_" + ("ok" if good else "REVISAR"))
    if prev_dat_last is not None:                                                  # promedios de ventana
        d = dict(T=dat[0, 1] - prev_dat_last[1], Ep=dat[0, 2] - prev_dat_last[2], Lx=dat[0, 3] - prev_dat_last[3], rho=dat[0, 4] - prev_dat_last[4])
        good = all(abs(d[k]) <= TOL[k] for k in TOL)
        print(f"   borde de VENTANAS (ultima fila prev -> primera de esta): dT={d['T']:+.2f} K  dLx={d['Lx']:+.4f} A  drho={d['rho']:+.5f}  dEp={1000*d['Ep']:+.2f} meV/at  "
              f"(tol {TOL['T']} K / {TOL['Lx']} / {TOL['rho']} / {1000*TOL['Ep']:.1f} meV)  {'OK' if good else 'REVISAR'}")
        ok_all &= good; status.append("ventanas_" + ("ok" if good else "REVISAR"))
    # stderr y log
    errs = sorted(glob.glob(f"{D}/slurm_seg{seg}.*.err")); bad = []
    for e in errs:
        for l in open(e, errors="ignore"):
            if not NOISE.search(l) and not re.match(r"^\s*\d+\)\s", l) and "unload" not in l and "following modules" not in l.lower(): bad.append(l.strip())
    lg = open(lp, errors="ignore").read(); flags = [k for k in ("ERROR", "lost atoms", "Fix halt condition") if k in lg]
    print(f"   stderr limpio: {'SI' if not bad else 'NO -> ' + ' | '.join(bad[:3])}   |  halt/ERROR/lost atoms en el log: {'ninguno' if not flags else flags}")
    ok_all &= (not bad and not flags); status.append("stderr_" + ("ok" if not bad else "SUCIO")); status.append("log_" + ("ok" if not flags else "PROBLEMA"))
    if os.path.exists(CHECK):
        o = subprocess.run([sys.executable, CHECK, dp], capture_output=True, text=True).stdout; line = [l for l in o.splitlines() if l.startswith(("OK", "ALERTA", "FLAG", "SE")) or "senal" in l.lower()]
        print("   cristalizacion:", (line[-1] if line else o.strip().splitlines()[-1] if o.strip() else "(sin salida)")); status.append("cryst:" + ("OK" if line and line[-1].startswith("OK") else "REVISAR"))
    rows.append([seg] + status); prev_end = end; prev_dat_last = dat[-1]
# ---------------- figura y alarma
def load_dat_all():
    xs = [read_dat(datfile(s)) for s in a.segs if datfile(s)]; return np.vstack(xs) if xs else None
def load_csv(p):
    r = [l.strip().split(",") for l in open(p)][1:]; r = np.array(r, float); return r[:, 1], r[:, 2]              # T, Ep
def binT(T, E, w=10.0):
    e = np.arange(np.floor(T.min() / w) * w, np.ceil(T.max() / w) * w + w, w); o = []
    for lo, hi in zip(e[:-1], e[1:]):
        m = (T >= lo) & (T < hi)
        if m.any(): o.append((T[m].mean(), E[m].mean()))
    o = np.array(o); return o[:, 0], o[:, 1]
br = load_dat_all()
if br is not None:
    ref = {"2e10": f"{SA}/2e10/thermo_2e10.csv", "1e11": f"{SA}/1e11/thermo_1e11.csv", "1e12": f"{SA}/1e12/thermo_1e12.csv", "qgrad": f"{SA}/qgrad/thermo/thermo_qgrad.csv"}
    B = {k: binT(*load_csv(p)) for k, p in ref.items()}; Tb, Eb = binT(br[:, 1], br[:, 2])
    COL = {"2e10": "#2a78d6", "1e11": "#eb6834", "1e12": "#1baf7a", "qgrad": "#7b3fbf"}; T2, E2 = B["2e10"]
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    for k in B: ax[0].plot(*B[k], ".", ms=2.5, color=COL[k], label=k)
    ax[0].plot(Tb, Eb, "-", color="k", lw=1.2, label=f"rama {P}"); ax[0].set(xlim=(280, 800), ylim=(-4.85, -4.72), xlabel="T (K)", ylabel="Ep (eV/at)", title="Ep(T)"); ax[0].legend(fontsize=8)
    off = lambda T, E: 1000 * (E - np.interp(T, T2, E2))
    for k in ("1e11", "1e12", "qgrad"): ax[1].plot(B[k][0], off(*B[k]), color=COL[k], label=f"{k} - 2e10")
    ax[1].plot(Tb, off(Tb, Eb), "k", lw=1.4, label=f"rama {P} - 2e10"); ax[1].axhline(0, color="grey", lw=.6); ax[1].set(xlim=(280, 800), xlabel="T (K)", ylabel="Ep - Ep(2e10) (meV/at)", title="brecha respecto de 2e10"); ax[1].legend(fontsize=8)
    M_G = 0.143; wins = [(lo, lo + 50) for lo in range(300, 800, 50)]; sl = {}
    for lo, hi in wins:
        m = (Tb >= lo) & (Tb <= hi)
        if m.sum() >= 4: sl[lo] = 1000 * np.polyfit(Tb[m], Eb[m], 1)[0]
    s2 = {lo: 1000 * np.polyfit(T2[(T2 >= lo) & (T2 <= lo + 50)], E2[(T2 >= lo) & (T2 <= lo + 50)], 1)[0] for lo, _ in wins}
    ax[2].plot([l + 25 for l in s2], list(s2.values()), "o-", color=COL["2e10"], ms=4, label="2e10"); ax[2].plot([l + 25 for l in sl], list(sl.values()), "s-", color="k", ms=4, label=f"rama {P}")
    ax[2].axhspan(M_G * .95, M_G * 1.05, color="grey", alpha=.25, label="vidrio 2e10 +-5%"); ax[2].set(xlabel="T (K), centro de ventana de 50 K", ylabel="dEp/dT (meV/at/K)", title="pendiente local"); ax[2].legend(fontsize=8)
    for x in ax: x.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, f"monitor_Ep_{P}.png"), dpi=140)
    inband = [lo for lo in sorted(sl, reverse=True) if abs(sl[lo] / M_G - 1) <= 0.05]; ks = sorted(sl, reverse=True)
    print("\n== Ep(T) y alarma ==\npendientes locales de la rama (meV/at/K, ventanas de 50 K, de mayor a menor T):", {f"{lo}-{lo+50}": round(sl[lo], 3) for lo in ks})
    print(f"pendiente vitrea de referencia (2e10, 300-500 K): {M_G} +-5% -> [{M_G*.95:.4f}, {M_G*1.05:.4f}]")
    alarm = None
    for i in range(len(ks) - 1):                                                   # dos ventanas consecutivas (de alta a baja T) en la banda
        l1, l2 = ks[i], ks[i + 1]
        if abs(sl[l1] / M_G - 1) <= .05 and abs(sl[l2] / M_G - 1) <= .05:
            mm = [(Tb >= lo) & (Tb <= lo + 50) for lo in (l1, l2)]; o = np.concatenate([off(Tb[m], Eb[m]) for m in mm])
            if o.max() - o.min() < 0.3: alarm = (l2, l1, o.max() - o.min()); break
    print("ALARMA 'rama lineal y paralela al vidrio de 2e10':", f"SI en {alarm[0]}-{alarm[1]+50} K (deriva de la brecha {alarm[2]:.2f} meV/at). Puede evaluarse cortar antes; la decision es tuya." if alarm else "no (todavia)")
    i0 = int(np.argmax(Tb >= 310)) if (Tb >= 310).any() else 0          # se excluye el bin de 300 K: el CSV de 2e10 incluye su relajacion
    print(f"brecha actual Ep - Ep(2e10) en el T mas bajo alcanzado por la rama ({Tb[i0]:.0f} K, bin de 300 K excluido): {off(Tb[i0:i0+1], Eb[i0:i0+1])[0]:+.2f} meV/at")
print("\nRESUMEN:", "todo OK" if ok_all else "HAY ITEMS A REVISAR (ver arriba)")
