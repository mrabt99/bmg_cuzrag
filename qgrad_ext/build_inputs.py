#!/usr/bin/env python3
"""Genera los inputs y sbatch de la rama qgrad_ext A PARTIR de los de qgrad (rate_gradiente_gammamin/), sin editar los originales.
Cada reemplazo se verifica con assert (si un patron no aparece, el script falla en vez de generar un input a medias).
Plan A: seg09 700->600 y seg10 600->500 a 1.5e9 K/s (copias de seg08 de qgrad); seg11 500->400 y seg12 400->300 a 1e11 K/s (copias de seg11/12 de qgrad);
relajacion 300 K/10 ns (copia de in.relax300K_qgrad). dt, ensemble (NPT iso 0 0, Tdamp 0.1, Pdamp 10), potencial, fix halt y frecuencias de thermo/ave-time/dump
NO se tocan: solo cambian T de la rampa, nombres de archivos (rate=gradext) y el restart de entrada.
Uso: python3 build_inputs.py"""
import hashlib, shutil
from pathlib import Path
Q = Path("/scratch/negishi/aaabt/bmg_cuzrag/rate_gradiente_gammamin"); R = Path("/scratch/negishi/aaabt/bmg_cuzrag/structural_analysis/qgrad/relax300K")
HERE = Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rep(s, old, new, n=None):
    assert old in s, f"patron no encontrado: {old[:70]!r}"
    return s.replace(old, new) if n is None else s.replace(old, new, n)
NOTICE = ("#  RAMA qgrad_ext (plan A): copia derivada de {src} (qgrad). Extiende/ajusta el perfil por debajo de 700 K.\n"
          "#  Solo cambian: rampa de T, nombres de archivos (rate=gradext) y restart de entrada. dt, ensemble, Tdamp/Pdamp, potencial,\n"
          "#  fix halt y frecuencias de thermo/ave-time/dump son IDENTICOS a qgrad. Generado por build_inputs.py.\n")
def notice(s, src): i = s.index("\n") + 1; return s[:i] + NOTICE.format(src=src) + s[i:]

# --- restart de partida (700 K, fin de la ventana lenta de qgrad) y potencial: copias con sha256 verificado
for src, dst in ((Q / "restart.quench_N65536_gradmin.seg08", HERE / "restart.start_700K_from_qgrad_seg8"), (Q / "ZrCuAg.lammps.eam", HERE / "ZrCuAg.lammps.eam")):
    if not dst.exists(): shutil.copy2(src, dst)
    assert sha(src) == sha(dst), f"sha256 distinto: {src}"; print(f"sha256 OK  {dst.name}  {sha(dst)[:16]}...")

# --- tramos lentos (base: in.quench_grad_seg08.lmp). seg10 se parte en 10a (600->550) y 10 (550->500) para tener restart en 550 K (opcion de cortar ahi).
# Mismas frecuencias de thermo/ave-time/dump que seg08 de qgrad (misma ventana de promedio); en los tramos de 33.3 M pasos quedan ~500 filas de log.dat y 6 frames.
base = (Q / "in.quench_grad_seg08.lmp").read_text()
SLOW = [("09", 700, 600, "restart.start_700K_from_qgrad_seg8", 66666667), ("10a", 600, 550, "restart.quench_N65536_gradext.seg09", 33333333),
        ("10", 550, 500, "restart.quench_N65536_gradext.seg10a", 33333334)]
for name, Ti, Tf, rin, steps in SLOW:
    s = base
    s = rep(s, "restart.quench_N65536_gradmin.seg07", rin)                         # comentario 'Entrada:' y read_restart
    s = rep(s, "variable        rate  index gradmin", "variable        rate  index gradext")
    s = rep(s, "seg08", f"seg{name}")                                              # log.dat / dump / write_restart / comentarios de salida
    s = rep(s, "tramo 8/12 : 800 K -> 700 K", f"tramo {name} (RAMA ext) : {Ti} K -> {Tf} K")
    s = rep(s, "velocity all create 800.0", f"velocity all create {Ti}.0")
    s = rep(s, '"### tramo 8/12 (gradmin, 1.500e+09 K/s) :  800 K -> 700 K"', f'"### tramo {name} ext (gradext, 1.500e+09 K/s) :  {Ti} K -> {Tf} K"')
    s = rep(s, "temp 800.0 700.0 0.1 iso 0.0 0.0", f"temp {Ti}.0 {Tf}.0 0.1 iso 0.0 0.0")
    s = rep(s, '"### tramo 8/12 terminado', f'"### tramo {name} ext terminado')
    if steps != 66666667:
        s = rep(s, "100K / tasa = 66666667 pasos", f"{Ti - Tf}K / tasa = {steps} pasos"); s = rep(s, "(66.667 ns)", f"({steps / 1e6:.3f} ns)")
        s = rep(s, "run              66666667                        # 66.667 ns", f"run              {steps}                        # {steps / 1e6:.3f} ns")
    (HERE / f"in.quench_ext_seg{name}.lmp").write_text(notice(s, "in.quench_grad_seg08.lmp"))
# --- tramos rapidos (base: seg11 y seg12 de qgrad)
for n in (11, 12):
    s = (Q / f"in.quench_grad_seg{n}.lmp").read_text()
    assert s.count("gradmin") >= 3
    s = s.replace("gradmin", "gradext")
    s = rep(s, f"tramo {n}/12 terminado", f"tramo {n} ext terminado"); s = rep(s, f"tramo {n}/12 (gradext", f"tramo {n} ext (gradext")
    (HERE / f"in.quench_ext_seg{n}.lmp").write_text(notice(s, f"in.quench_grad_seg{n}.lmp"))
# --- relajacion 300 K (base: in.relax300K_qgrad de qgrad/relax300K)
s = (R / "in.relax300K_qgrad").read_text()
s = rep(s, "variable        pot   index /scratch/negishi/aaabt/bmg_cuzrag/rate_gradiente_gammamin/ZrCuAg.lammps.eam", "variable        pot   index ZrCuAg.lammps.eam")
s = rep(s, "variable        rate  index qgrad", "variable        rate  index gradext")
s = rep(s, "read_restart    /scratch/negishi/aaabt/bmg_cuzrag/rate_gradiente_gammamin/restart.quench_N65536_gradmin.seg12", "read_restart    restart.quench_N65536_gradext.seg12")
s = rep(s, "#  Arranca de restart.quench_N65536_gradmin.seg12 (301 K, sin relajar).", "#  RAMA qgrad_ext: arranca de restart.quench_N65536_gradext.seg12 (~300 K, sin relajar). Copia de in.relax300K_qgrad.")
s = rep(s, "#  NO modifica nada de rate_gradiente_gammamin/ (solo lee el restart);", "#  Todo queda en qgrad_ext/.")
s = rep(s, "#  todo lo generado queda en structural_analysis/qgrad/relax300K/.", "#  (mismo NPT 300 K / 0 bar / 10 ns / Pdamp 10 / Tdamp 0.1 / dt 1 fs que el 17r de qgrad y las otras tasas)")
(HERE / "in.relax300K_ext").write_text(s)
# --- sbatch (bases: slurm_seg08 / seg11 / seg12 de qgrad y slurm_relax de qgrad/relax300K)
def sbatch(src, n, inp, name, walltime=None, extra=()):
    t = src.read_text()
    return t
for name, Ti, Tf, steps in (("09", 700, 600, 66666667), ("10a", 600, 550, 33333333), ("10", 550, 500, 33333334)):
    t = (Q / "slurm_seg08.sbatch").read_text()
    t = rep(t, "qgrad_seg08", f"qgradext_seg{name}"); t = rep(t, "slurm_seg08", f"slurm_seg{name}"); t = rep(t, "in.quench_grad_seg08.lmp", f"in.quench_ext_seg{name}.lmp"); t = rep(t, "log.seg08", f"log.seg{name}")
    t = rep(t, "tramo 8/12 (800->700K, 1.500e+09 K/s)", f"RAMA ext, tramo {name} ({Ti}->{Tf} K, 1.500e+09 K/s)")
    if steps != 66666667:     # tramo de la mitad de largo: 33.3 M pasos =~ 30.9 h a 300 pasos/s -> 2.2x = 2-19:54:00
        t = rep(t, "66666667 pasos =~ 61.7h", f"{steps} pasos =~ 30.9h"); t = rep(t, "5-16:19:00 de margen", "2-19:54:00 de margen"); t = rep(t, "#SBATCH --time=5-16:19:00", "#SBATCH --time=2-19:54:00")
    (HERE / f"slurm_seg{name}.sbatch").write_text(t)
for n in (11, 12):
    t = (Q / f"slurm_seg{n}.sbatch").read_text()
    t = rep(t, f"qgrad_seg{n}", f"qgradext_seg{n}"); t = rep(t, f"in.quench_grad_seg{n}.lmp", f"in.quench_ext_seg{n}.lmp")
    (HERE / f"slurm_seg{n}.sbatch").write_text(t)
t = (R / "slurm_relax.sbatch").read_text()
t = rep(t, "relax300K_qgrad", "relax300K_ext")     # cubre tambien "in.relax300K_qgrad" -> "in.relax300K_ext"
(HERE / "slurm_relax.sbatch").write_text(t)
print("inputs y sbatch generados:", sorted(p.name for p in HERE.glob("in.*")) + sorted(p.name for p in HERE.glob("slurm_*.sbatch")))
