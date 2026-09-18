#!/usr/bin/env python3
import sys, math

path = sys.argv[1] if len(sys.argv) > 1 else "ZrCuAg.lammps.eam"
lines = open(path).readlines()
st = {"idx": 3, "buf": []}

def header_line():
    st["buf"] = []
    l = lines[st["idx"]]; st["idx"] += 1
    return l.split('#')[0].split()

def take(n):
    buf = st["buf"]
    while len(buf) < n:
        l = lines[st["idx"]]; st["idx"] += 1
        buf.extend(l.split('#')[0].split())
    out = [float(v) for v in buf[:n]]
    del buf[:n]
    return out

parts = header_line()
nelem = int(parts[0]); elements = parts[1:1+nelem]
p5 = header_line()
nrho, drho = int(p5[0]), float(p5[1])
nr, dr, cutoff = int(p5[2]), float(p5[3]), float(p5[4])

Z = {}; F = {}
for e in elements:
    h = header_line()
    Z[e] = int(float(h[0]))
    F[e] = take(nrho)
    take(nr)                      # rho(r), no lo necesitamos

phi = {}
for i in range(nelem):
    for j in range(i+1):
        phi[(elements[i], elements[j])] = take(nr)

def zbl(z1, z2, rr):
    a = 0.46850/(z1**0.23 + z2**0.23); x = rr/a
    s = (0.18175*math.exp(-3.19980*x) + 0.50986*math.exp(-0.94229*x)
         + 0.28022*math.exp(-0.40290*x) + 0.02817*math.exp(-0.20162*x))
    return 14.399645*z1*z2/rr*s

print("elementos %s  Z=%s" % (elements, [Z[e] for e in elements]))
print("nr=%d dr=%g  nrho=%d drho=%g  cutoff=%g\n" % (nr, dr, nrho, drho, cutoff))

print("--- phi(r) en eV: EAM vs ZBL ---")
for (a, b), rphi in phi.items():
    print("%s-%s" % (a, b))
    for x in [0.4, 0.6, 0.8, 1.0, 1.5]:
        k = int(round(x/dr))
        if k == 0 or k >= nr:
            continue
        e_val = rphi[k]/(k*dr)
        z_val = zbl(Z[a], Z[b], x)
        print("    r=%.1f: EAM %10.1f   ZBL %9.1f   x%.3f" % (x, e_val, z_val, e_val/z_val))

print("\n--- F(rho) en la mitad alta de la tabla (debe crecer siempre) ---")
for e in elements:
    Fe = F[e]; half = nrho//2
    bad = [k for k in range(half, nrho-1) if Fe[k+1] - Fe[k] < 0]
    if bad:
        print("  %s: PROBLEMA, dF/drho < 0 en rho ~ %.3f a %.3f"
              % (e, bad[0]*drho, bad[-1]*drho))
    else:
        print("  %s: ok, monotona creciente hasta rho_max=%.3f" % (e, (nrho-1)*drho))
