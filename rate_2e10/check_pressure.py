#!/usr/bin/env python3
"""Monitor de presion/temperatura para los log.dat de fix ave/time del quench 1e10 K/s.

Detecta columnas por nombre EXACTO de header (no substring): en v3/sol67 se
encontro que buscar "P(" con `in` da falso positivo porque "P(" es substring
de "Ep(". Aca se parte el header en tokens y se compara por igualdad.

Uso:
    python3 check_pressure.py quench_N65536_1e10.log.dat.seg1 [seg2 ...]
    python3 check_pressure.py --fracture 2.3 quench_N65536_1e10.log.dat.seg*
"""
import argparse
import glob
import sys

BAR_PER_GPA = 1.0e4  # unidades metal de LAMMPS: presion en bar


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
                # el titulo real (title2) es la ultima linea de comentario
                # antes de los datos; nos quedamos con la mas reciente.
                # OJO: fix ave/time siempre escribe el numero de paso como
                # primera columna de datos, pero title2 en estos inputs
                # solo etiqueta las columnas v_* (no incluye "step") ->
                # hay que anteponerlo a mano o los indices quedan corridos
                # uno a la izquierda (bug real: sin esto, idx_T apuntaba a
                # la columna de step y idx_P a la de rho).
                header_tokens = ["step"] + line.lstrip("#").split()
                continue
            rows.append(line.split())
    return header_tokens, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", help="archivos log.dat (o patrones glob)")
    ap.add_argument("--fracture", type=float, default=2.3,
                     help="umbral de alerta en GPa (default 2.3 GPa, sigma_y de "
                          "Mendez & Jaroszewicz)")
    args = ap.parse_args()

    paths = []
    for pat in args.files:
        matches = sorted(glob.glob(pat))
        paths.extend(matches if matches else [pat])

    worst_gpa = None
    worst_loc = None
    any_over = False

    for path in paths:
        header_tokens, rows = parse_file(path)
        if not header_tokens:
            print(f"[{path}] sin header reconocible, salteando", file=sys.stderr)
            continue

        # title2 en los inputs de este proyecto es:
        # "# T(K)  Ep(eV/at)  Lx(A)  rho(g/cc)  P(bar)"
        idx_T = find_column(header_tokens, "T(K)")
        idx_P = find_column(header_tokens, "P(bar)")
        if idx_T is None or idx_P is None:
            print(f"[{path}] no encontre columnas T(K)/P(bar) exactas en "
                  f"header={header_tokens}, salteando", file=sys.stderr)
            continue

        for row in rows:
            if len(row) <= max(idx_T, idx_P):
                continue
            try:
                step = row[0]
                T = float(row[idx_T])
                P_bar = float(row[idx_P])
            except ValueError:
                continue
            P_gpa = P_bar / BAR_PER_GPA
            if worst_gpa is None or abs(P_gpa) > abs(worst_gpa):
                worst_gpa = P_gpa
                worst_loc = (path, step, T)
            if abs(P_gpa) > args.fracture:
                any_over = True
                print(f"[{path}] step={step} T={T:.1f}K  |P|={abs(P_gpa):.3f} GPa "
                      f"> umbral {args.fracture} GPa")

    if worst_loc:
        path, step, T = worst_loc
        pct = 100.0 * abs(worst_gpa) / args.fracture
        print(f"\nPico de |P| = {abs(worst_gpa):.3f} GPa en {path} "
              f"(step={step}, T={T:.1f} K) -> {pct:.1f}% del umbral "
              f"({args.fracture} GPa)")

    if any_over:
        print("\nALERTA: hubo puntos por encima del umbral de fractura.")
        sys.exit(1)
    else:
        print("\nOK: ningun punto supero el umbral de fractura.")


if __name__ == "__main__":
    main()
