"""Fase 5: tablas finales de transporte (LaTeX + CSV).
  tabla_kappa.tex      kappa(T) por tasa (T medida; error = sem de bloques (+) sistematico de ventana)
  tabla_eta_D.tex      eta, D y Stokes-Einstein en el liquido
  tabla_corridas.csv/.tex  parametros de cada corrida, leidos de los archivos de cada carpeta cal_*
Fuentes: ../kappa_tasas.csv, ../eta_T.csv, ../eta_tasas.csv, ../eta_largo.log (900 y 800 K, 2e10).
Semillas: convencion de los lanzadores (r1-r3: 1001/2002/3003; r4-r6: 4004/5005/6006), no estan en los logs.
Uso (desde transporte/fase5/): python tablas.py"""
import csv, glob, os, re, numpy as np
rd=lambda f: list(csv.DictReader(open(f)))
KT=rd("../kappa_tasas.csv"); ET=rd("../eta_T.csv"); EA=rd("../eta_tasas.csv"); LOG=open("../eta_largo.log").read()
RATES=["2e10","1e11","1e12"]; HDR={"2e10":r"$2\times10^{10}$","1e11":r"$10^{11}$","1e12":r"$10^{12}$"}
def fmt(v,e,n=3): return f"${v:.{n}f}\\pm{e:.{n}f}$"
# ---------- tabla kappa ----------
L=[r"\begin{table}[htbp]\centering",
   r"\caption{Conductividad térmica de red $\kappa$ (W\,m$^{-1}$\,K$^{-1}$) por Green--Kubo, integral por trapecio promediada en la ventana 0.5--2~ps. "
   r"Tres réplicas de 1~ns por punto (30 bloques de 100~ps). El error combina el error estándar entre bloques y la dispersión entre ventanas. "
   r"$^\dagger$: $\kappa(t)$ no pasa el test de meseta ($|$pendiente$|\geq 2\sigma$). Entre 1300 y 1600~K solo se midió la tasa $2\times10^{10}$~K/s (error solo entre bloques).}",
   r"\label{tab:kappa}",r"\begin{tabular}{rccc}",r"\hline",r"$T$ (K) & "+" & ".join(HDR[r]+" K/s" for r in RATES)+r" \\",r"\hline"]
for T in sorted({int(x["T_K"]) for x in KT}):
    c=[]
    for r in RATES:
        x=[y for y in KT if int(y["T_K"])==T and y["tasa"]==r][0]
        c.append(fmt(float(x["kappa"]),np.hypot(float(x["sem"]),float(x["sistematico"])))+(r"$^\dagger$" if x["meseta"]=="no" else ""))
    L.append(f"{T} & "+" & ".join(c)+r" \\")
for x in ET:
    if int(x["T_K"])>=1300: L.append(f"{x['T_K']} & "+fmt(float(x["kappa"]),float(x["kappa_sem"]))+r" & -- & -- \\")
L+=[r"\hline",r"\end{tabular}",r"\end{table}"]; open("tabla_kappa.tex","w").write("\n".join(L)+"\n")
# ---------- tabla eta, D ----------
def bl(T): return LOG.split(f"===== {T} K")[1].split("=====")[0]
D900=[float(v) for v in re.search(r"D\(Zr,Cu,Ag\) = ([\d.]+) ([\d.]+) ([\d.]+)",bl(900)).groups()]
E900=[float(v) for v in re.search(r"ventana \[20,60\] ps: eta = ([\d.]+) \+- ([\d.]+)",bl(900)).groups()]
D800=[float(v) for v in re.search(r"D\(Zr,Cu,Ag\) = ([\d.]+) ([\d.]+) ([\d.]+)",bl(800)).groups()]
E800=[float(v) for v in re.search(r"t=199:([\d.]+)\+-([\d.]+)",bl(800)).groups()]
rows={800:(f"$>{E800[0]:.0f}$",*D800),900:(fmt(*E900,1),*D900)}
for x in ET: rows[int(x["T_K"])]=(fmt(float(x["eta_mPas"]),float(x["sem"]),(1 if float(x["eta_mPas"])>=10 else 2)),float(x["D_Zr"]),float(x["D_Cu"]),float(x["D_Ag"]))
def eo(r,T):
    q=[y for y in EA if y["tasa"]==r and int(y["T_K"])==T]
    if not q or T<1000: return "--"
    v=float(q[0]["eta_Tnom"]); e=float(q[0]["sem"]); return fmt(v,e,(1 if v>=10 else 2))
L=[r"\begin{table}[htbp]\centering",
   r"\caption{Viscosidad de corte $\eta$ (mPa\,s) y coeficientes de difusión $D$ ($10^{-9}$\,m$^2$\,s$^{-1}$) en el líquido. "
   r"$2\times10^{10}$~K/s: 3 réplicas de 1~ns (6 a 1000~K); a 900~K, 3 réplicas de 3~ns y ventana 20--60~ps; a 800~K, 3 réplicas de 5~ns, donde $\eta(t)$ no alcanza meseta y se da la cota $\eta(200~\mathrm{ps})$. "
   r"$10^{11}$ y $10^{12}$~K/s: 3 réplicas de 1~ns, valores llevados a la temperatura nominal. A 900~K, con 1~ns por réplica, $\eta(20~\mathrm{ps})$ vale 37.1, 37.5, 38.2 en las tres tasas (subestima en $\sim$12\,\% el valor con 3~ns). $D$ por especie a $2\times10^{10}$~K/s, de la pendiente del desplazamiento cuadrático medio.}",
   r"\label{tab:eta}",r"\begin{tabular}{rccccccc}",r"\hline",
   r"& \multicolumn{3}{c}{$\eta$} & \multicolumn{3}{c}{$D$ ($2\times10^{10}$~K/s)} \\",
   r"$T$ (K) & "+" & ".join(HDR[r] for r in RATES)+r" & Zr & Cu & Ag \\",r"\hline"]
for T in sorted(rows):
    e,dz,dc,da=rows[T]; L.append(f"{T} & {e} & {eo('1e11',T)} & {eo('1e12',T)} & {dz:.3f} & {dc:.3f} & {da:.3f}"+r" \\")
L+=[r"\hline",r"\end{tabular}",r"\end{table}"]; open("tabla_eta_D.tex","w").write("\n".join(L)+"\n")
# ---------- tabla de corridas (leida de los archivos) ----------
out=[]
for d in sorted(glob.glob("../cal_*K_*"),key=lambda s:(("N8192" in s),int(re.search(r"cal_(\d+)K",s).group(1)),s)):
    name=os.path.basename(d); T=int(re.search(r"cal_(\d+)K",name).group(1)); rate=name.split("_")[2]
    npt=open(f"{d}/run_npt.sbatch").read() if os.path.exists(f"{d}/run_npt.sbatch") else ""
    m=re.search(r"-var rst (\S+)",npt); rst=os.path.basename(m.group(1)) if m else ""
    if not rst and os.path.exists(f"{d}/in.npt"):
        m=re.search(r"variable\s+rst\s+index\s+(\S+)",open(f"{d}/in.npt").read()); rst=os.path.basename(m.group(1)) if m else "?"
    m=re.search(r"-var nsteps (\d+)",npt); nnpt=int(m.group(1)) if m else int(re.search(r"variable\s+nsteps\s+index\s+(\d+)",open(f"{d}/in.npt").read()).group(1))
    nrep=len(glob.glob(f"{d}/cal_r*_hflux.dat")); nv=open(f"{d}/run_nve.sbatch").read(); nprod=int(re.search(r"-var nprod (\d+)",nv).group(1)) if "-var nprod" in nv else 0
    rows_=[(name,T,rate,"8192" if "N8192" in name else "65536",rst,nnpt/1000,f"r1-r{nrep}",nrep,nprod/1e6,"1001,2002,...",2,5,"0.5-2 ps","[t1,3 t1]" if T>=1000 else ("--" if T<900 else "[t1,3 t1] o eta(20 ps)"))]
    nl=sorted({re.search(r"long_r(\d+)_p",f).group(1) for f in glob.glob(f"{d}/long_r*_p*_hflux.dat")})
    if nl:
        npart=len(glob.glob(f"{d}/long_r{nl[0]}_p*_hflux.dat"))
        rows_.append((name+" (largas)",T,rate,"65536",rst,nnpt/1000,f"r{nl[0]}-r{nl[-1]}",len(nl),float(npart),"4004,5005,6006",2,5,"0.5-2 ps","[20,60] ps" if T==900 else "sin meseta (cota)"))
    out+=rows_
H=["carpeta","T_K","tasa","N","restart","NPT_ps","replicas","n_rep","NVE_ns_por_rep","semillas","muestreo_P_fs","muestreo_J_fs","ventana_kappa","ventana_eta"]
with open("tabla_corridas.csv","w") as f:
    f.write(",".join(H)+"\n"); [f.write(",".join(str(v).replace(",",";") for v in r)+"\n") for r in out]
L=[r"\begingroup\footnotesize",r"\begin{longtable}{llrlrcrl}",
   r"\caption{Corridas de Green--Kubo. Todas: $\Delta t=1$~fs, NPT con $\tau_T=0.1$~ps y $\tau_P=10$~ps, NVT de 20~ps por réplica antes del NVE, presión muestreada cada 2~fs y flujo de calor cada 5~fs. "
   r"Las réplicas 1--3 usan las semillas 1001, 2002 y 3003; las 4--6, 4004, 5005 y 6006.}\label{tab:corridas}\\",
   r"\hline",r"$T$ (K) & tasa & $N$ & restart & NPT (ps) & réplicas & NVE (ns) & ventana $\eta$ \\",r"\hline\endhead"]
for r in out:
    L.append(f"{r[1]} & {r[2]} & {r[3]} & \\texttt{{{r[4].replace('restart.','').replace('_N65536','').replace('_',chr(92)+'_')}}} & {r[5]:.0f} & {r[7]} & {r[8]:.0f} & {r[13].replace('t1',chr(36)+'t_1'+chr(36))}"+r" \\")
L+=[r"\hline",r"\end{longtable}",r"\endgroup"]; open("tabla_corridas.tex","w").write("\n".join(L)+"\n")
print(open("tabla_kappa.tex").read()); print(open("tabla_eta_D.tex").read()); print(f"tabla_corridas: {len(out)} filas"); print(open("tabla_corridas.csv").read())
