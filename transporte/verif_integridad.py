"""Chequea integridad de las salidas de fix ave/time (cal_r*_{hflux,press,thermo,msd}.dat) en los dirs dados:
bytes nulos, lineas con numero de columnas incorrecto, pasos de tiempo no monotonicos o con salto distinto al esperado,
encabezados repetidos en el medio. Uso: python verif_integridad.py dir1 dir2 ..."""
import sys, glob, numpy as np
for d in sys.argv[1:]:
    for f in sorted(glob.glob(f"{d}/cal_r*_*.dat")):
        raw=open(f,"rb").read(); nul=raw.count(b"\x00")
        txt=raw.replace(b"\x00",b"").decode(errors="replace").splitlines()
        hdr=[i for i,l in enumerate(txt) if l.startswith("#")]
        rows=[l.split() for l in txt if l.strip() and not l.startswith("#")]
        nc=max(set(len(r) for r in rows),key=[len(r) for r in rows].count)
        bad=sum(len(r)!=nc for r in rows)
        ts=np.array([int(float(r[0])) for r in rows if len(r)==nc and r[0].lstrip("-").replace(".","").isdigit()])
        dts=np.diff(ts); step=np.median(dts) if len(dts) else 0
        nonmono=int((dts!=step).sum())
        ok = nul==0 and bad==0 and nonmono==0 and hdr==[0,1]
        print(f"{'OK ' if ok else 'MAL'} {f}: filas={len(rows)} ts {ts[0] if len(ts) else '-'}..{ts[-1] if len(ts) else '-'} paso={step:.0f}"
              f" nulos={nul} cols_mal={bad} saltos_raros={nonmono} encabezados_en={hdr[:4]}")
