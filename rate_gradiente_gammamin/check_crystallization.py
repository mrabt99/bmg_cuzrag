#!/usr/bin/env python3
"""Monitor de cristalizacion para la ventana graduada 1500->700K (gamma_min).

Ryltsev et al. reportan cristalizacion parcial (Cu2Zr) cerca de gamma_min
(1.5e9 K/s). Esta corrida es exploratoria y baja de 2e10 a 1.5e9 K/s en
8 tramos -- el chequeo mas importante NO es la explosion (eso ya lo cubre
fix halt), sino una caida ABRUPTA de energia potencial y/o un salto de
densidad por encima de la tendencia suave esperada del enfriamiento.

Metodo: baseline por mediana movil (ventana configurable) sobre PEat y
rho, residuo = valor - baseline, umbral en sigma robusto (MAD*1.4826)
sobre esos residuos. Flag cuando PEat cae Y rho sube mas de lo esperado
cerca del mismo punto (senal conjunta = mas confiable que una sola serie).

Uso:
    python3 check_crystallization.py quench_N65536_gradmin.log.dat.seg0*
    python3 check_crystallization.py --sigma 6 --window 15 log.dat.seg*
"""
import argparse
import glob
import statistics
import sys


def find_column(header_tokens, name):
    for idx, tok in enumerate(header_tokens):
        if tok == name:
            return idx
    return None


def parse_file(path):
    header_tokens = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # mismo gotcha que check_pressure.py: title2 no incluye la
                # columna de step que fix ave/time siempre antepone.
                header_tokens = ["step"] + line.lstrip("#").split()
                continue
            rows.append(line.split())
    return header_tokens, rows


def robust_sigma(residuals):
    med = statistics.median(residuals)
    mad = statistics.median([abs(r - med) for r in residuals])
    return mad * 1.4826 if mad > 0 else None


def rolling_median_baseline(values, window):
    n = len(values)
    half = window // 2
    baseline = []
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        baseline.append(statistics.median(values[lo:hi]))
    return baseline


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", help="archivos log.dat.segNN, en orden")
    ap.add_argument("--sigma", type=float, default=8.0,
                     help="umbral en sigma robusto (default 8, ver memoria del "
                          "proyecto 1e11: sin discontinuidades >8sigma = sin senal)")
    ap.add_argument("--window", type=int, default=21,
                     help="ventana de la mediana movil (en puntos de salida, no pasos)")
    args = ap.parse_args()

    paths = []
    for pat in args.files:
        matches = sorted(glob.glob(pat))
        paths.extend(matches if matches else [pat])

    all_rows = []  # (path, step, T, PEat, Lx, Rho, Pbar)
    for path in paths:
        header_tokens, rows = parse_file(path)
        if not header_tokens:
            print(f"[{path}] sin header reconocible, salteando", file=sys.stderr)
            continue
        idx = {name: find_column(header_tokens, name)
               for name in ("T(K)", "Ep(eV/at)", "Lx(A)", "rho(g/cc)", "P(bar)")}
        if any(v is None for v in idx.values()):
            print(f"[{path}] faltan columnas exactas en header={header_tokens}, "
                  f"salteando", file=sys.stderr)
            continue
        for row in rows:
            if len(row) <= max(idx.values()):
                continue
            try:
                step = int(row[0])
                T = float(row[idx["T(K)"]])
                PEat = float(row[idx["Ep(eV/at)"]])
                Lx = float(row[idx["Lx(A)"]])
                Rho = float(row[idx["rho(g/cc)"]])
                Pbar = float(row[idx["P(bar)"]])
            except ValueError:
                continue
            all_rows.append((path, step, T, PEat, Lx, Rho, Pbar))

    if len(all_rows) < args.window * 2:
        print(f"Muy pocos puntos ({len(all_rows)}) para un baseline confiable "
              f"con ventana {args.window} -- ejecutar con mas segmentos o "
              f"bajar --window.", file=sys.stderr)
        sys.exit(2)

    PEat_series = [r[3] for r in all_rows]
    Rho_series = [r[5] for r in all_rows]

    PEat_base = rolling_median_baseline(PEat_series, args.window)
    Rho_base = rolling_median_baseline(Rho_series, args.window)

    PEat_resid = [v - b for v, b in zip(PEat_series, PEat_base)]
    Rho_resid = [v - b for v, b in zip(Rho_series, Rho_base)]

    sigma_PE = robust_sigma(PEat_resid)
    sigma_Rho = robust_sigma(Rho_resid)

    flags = []
    for i, row in enumerate(all_rows):
        path, step, T, PEat, Lx, Rho, Pbar = row
        z_pe = PEat_resid[i] / sigma_PE if sigma_PE else 0.0
        z_rho = Rho_resid[i] / sigma_Rho if sigma_Rho else 0.0
        # cristalizacion: PE cae por debajo de la tendencia (z_pe muy negativo)
        # Y rho sube por encima de la tendencia (z_rho muy positivo), cerca
        # del mismo punto -- exigimos la senal conjunta para evitar falsos
        # positivos por ruido termico normal en una sola serie.
        if z_pe < -args.sigma and z_rho > args.sigma:
            flags.append((path, step, T, PEat, Rho, z_pe, z_rho))

    print(f"Puntos analizados: {len(all_rows)}  (sigma_PE={sigma_PE:.3e} eV/at, "
          f"sigma_rho={sigma_Rho:.3e} g/cc, umbral={args.sigma} sigma)")

    if flags:
        print(f"\nALERTA: {len(flags)} punto(s) con caida conjunta de PE + "
              f"salto de densidad > {args.sigma} sigma (posible cristalizacion):")
        for path, step, T, PEat, Rho, z_pe, z_rho in flags:
            print(f"  [{path}] step={step} T={T:.1f}K  PEat={PEat:.5f}eV "
                  f"(z={z_pe:+.1f})  rho={Rho:.4f}g/cc (z={z_rho:+.1f})")
        print("\nNo sigas la cadena mas alla de este tramo sin avisar al usuario.")
        sys.exit(1)
    else:
        print("\nOK: sin senal conjunta de cristalizacion en este rango.")


if __name__ == "__main__":
    main()
