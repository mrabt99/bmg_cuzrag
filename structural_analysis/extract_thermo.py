#!/usr/bin/env python3
"""Etapa 1a: concatena la historia termica completa de una tasa de
enfriamiento (todos los segmentos, no solo los 17 snapshots) en una
tabla continua: step, T, Ep(eV/at), Lx, rho, P.

Dos formatos de entrada, porque las tasas NO comparten el mismo tipo
de log (ver discusion con el usuario, 2026-09-16):

  --style avetime   2e10 / 1e11: N archivos .log.dat de `fix ave/time`,
                    uno por segmento, en orden fisico. Cada archivo
                    tiene 2 lineas de comentario; la 2da (title2) NO
                    nombra la primera columna real (el paso que
                    ave/time antepone siempre) -> por eso los nombres
                    de columna estan hardcodeados aca (step, T, Ep,
                    Lx, rho, P), nunca leidos del title2.

  --style rawthermo 1e12: un unico log de texto crudo de LAMMPS (thermo
                    screen output, `Step Temp Press Volume Density
                    PotEng TotEng`) con un bloque por cada `run`. Acá
                    PotEng es TOTAL (units metal => thermo_modify norm
                    no por defecto), hay que dividir por N para el
                    Ep/atomo; y no hay columna Lx, se deriva de
                    Volume**(1/3) (valido porque todo el proyecto usa
                    NPT iso => caja cubica). Los bloques consecutivos
                    comparten el timestep de arranque (cada `run`
                    repite el thermo del ultimo paso del run anterior)
                    -> se deduplica por step.
"""
import argparse
import csv


def read_avetime_logs(log_files):
    """Concatena N archivos fix-ave/time (step T Ep Lx rho P), en orden."""
    rows = []
    for path in log_files:
        with open(path) as f:
            lines = [l for l in f if not l.startswith("#")]
        for line in lines:
            parts = line.split()
            if len(parts) != 6:
                continue
            step, T, Ep, Lx, rho, P = (float(x) for x in parts)
            rows.append((int(step), T, Ep, Lx, rho, P))
    return rows


def read_rawthermo_log(stdout_file, natoms):
    """Parsea un log crudo de LAMMPS con multiples bloques thermo_style
    custom `step temp press vol density pe etotal`, deduplicando el
    step de arranque compartido entre bloques consecutivos."""
    header = None
    blocks = []
    current = None
    with open(stdout_file) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("Step Temp Press Volume Density PotEng"):
                if current:
                    blocks.append(current)
                current = []
                continue
            if current is None:
                continue
            parts = stripped.split()
            if len(parts) != 7:
                # linea no-numerica (Loop time..., blank, etc.) -> fin de bloque
                if current:
                    blocks.append(current)
                current = None
                continue
            try:
                step, T, P, vol, rho, pot, tot = (float(x) for x in parts)
            except ValueError:
                if current:
                    blocks.append(current)
                current = None
                continue
            current.append((int(step), T, P, vol, rho, pot, tot))
    if current:
        blocks.append(current)

    rows = []
    last_step = None
    for block in blocks:
        for (step, T, P, vol, rho, pot, tot) in block:
            if step == last_step:
                continue  # boundary row repetida entre bloques consecutivos
            Lx = vol ** (1.0 / 3.0)
            Ep = pot / natoms
            rows.append((step, T, Ep, Lx, rho, P))
            last_step = step
    return rows


def write_csv(rows, out_csv):
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "T", "Ep", "Lx", "rho", "P"])
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {out_csv}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--style", required=True, choices=["avetime", "rawthermo"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--logs", nargs="+", help="avetime: lista ordenada de .log.dat")
    ap.add_argument("--stdout", help="rawthermo: archivo unico de log crudo")
    ap.add_argument("--natoms", type=int, default=65536)
    args = ap.parse_args()

    if args.style == "avetime":
        if not args.logs:
            ap.error("--style avetime requiere --logs")
        rows = read_avetime_logs(args.logs)
    else:
        if not args.stdout:
            ap.error("--style rawthermo requiere --stdout")
        rows = read_rawthermo_log(args.stdout, args.natoms)

    write_csv(rows, args.out)


if __name__ == "__main__":
    main()
