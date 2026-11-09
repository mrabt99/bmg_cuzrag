#!/usr/bin/env python3
"""Verifica la integridad de las salidas de TEXTO de un caso de traccion (curva, dump, log): sin bytes NUL, 200 lineas de curva con 11 columnas, 41 cuadros de dump, y los 5 restarts con el tamano esperado.
(Los .restart son binarios: tienen NUL legitimos, solo se comprueba el tamano.) Uso: python check_integrity.py <carpeta del caso>   o   from check_integrity import check(path) -> lista de problemas ([] = integro)
Motivo: el 2026-10-08 un incidente de E/S en /scratch dejo jobs 'COMPLETED' con bloques de NUL y cuadros de dump faltantes."""
import sys, os, glob, re
def check(p):
    p = os.path.realpath(p); prob = []
    cur = glob.glob(f"{p}/*.curve.dat"); dmp = glob.glob(f"{p}/*.dump"); logs = sorted(glob.glob(f"{p}/log.tens_*"))
    if not (cur and dmp and logs): return ["faltan curve/dump/log"]
    raw = open(cur[0], "rb").read(); ok = [l for l in raw.split(b"\n") if l and not l.startswith(b"#")]
    if raw.count(b"\x00") or any(len(l.split()) != 11 for l in ok): prob.append(f"curva con {raw.count(chr(0).encode())} NUL / lineas rotas")
    if len(ok) != 200: prob.append(f"curva {len(ok)} lineas (200)")
    data = open(dmp[0], "rb").read(); nf = len(re.findall(rb"ITEM: TIMESTEP\n(\d+)", data))
    if data.count(b"\x00"): prob.append(f"dump con {data.count(b'\x00')} NUL")
    if nf != 41: prob.append(f"dump {nf} cuadros (41)")
    lg = open(logs[-1], "rb").read()
    if lg.count(b"\x00"): prob.append(f"log con {lg.count(b'\x00')} NUL")
    if b"### TRACCION TERMINADA" not in lg: prob.append("log sin TRACCION TERMINADA")
    rs = [f for f in glob.glob(f"{p}/*.restart")]; big = [f for f in rs if not f.endswith("ref.restart")]
    if len(rs) != 5 or any(os.path.getsize(f) != 5770427 for f in big): prob.append("restarts incompletos o de tamano distinto")
    return prob
if __name__ == "__main__":
    pr = check(sys.argv[1]); print("INTEGRO" if not pr else "DANADO: " + "; ".join(pr)); sys.exit(1 if pr else 0)
