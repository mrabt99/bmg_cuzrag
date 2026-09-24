#!/usr/bin/env python3
"""Snapshot 3D del vidrio 300 K coloreado por clase de cluster (analogo a Fig 10, Ryltsev 2016).

Uso (desde el directorio que contiene ovito_in/):
    python ovito_snapshot.py                 # las 3 tasas
    python ovito_snapshot.py 2e10 1e12       # solo algunas
    python ovito_snapshot.py --no-render     # solo conteos y chequeos
    python ovito_snapshot.py --opengl        # renderer OpenGL (rapido, requiere GL/display)
    python ovito_snapshot.py --force         # sigue aunque falle el chequeo de CN
    python ovito_snapshot.py --slab 15       # corte de 15 A (normal z, centro de la caja), vista superior; archivos *_slab_*
                                             # (clasificacion sobre el sistema COMPLETO; el corte es solo visual)
    python ovito_snapshot.py --no-ico-voronoi # variante: SIN clase ico-Voronoi (esos atomos caen en "resto");
                                             # archivos con tag "noicovor_" (no pisa los de la version completa)
    python ovito_snapshot.py --bare          # render SIN texto -> <tasa>_bare.png (insumo de compose_snapshots.py)
    python ovito_snapshot.py --resto-transp 1 --radius-scale 0.65
                                             # variante "despejada": --resto-transp X pone la transparencia de "resto" en X
                                             # (X=1 => los atomos "resto" se eliminan del render); --radius-scale S multiplica el
                                             # radio de las clases coloreadas (no de "resto"). Solo visual: la clasificacion no cambia.
                                             # Archivos con tag "rt<X>_rs<S>_" (no pisan los de la version completa)
    python ovito_snapshot.py --center Cu     # solo dibuja los atomos CENTRALES de esa especie (Cu|Zr|Ag); los de otras especies
                                             # caen en "resto" (solo visual). counts_*.csv de esa variante: N por clase y n_total =
                                             # numero de atomos de la especie. Tag "<especie>_" (se combina con los demas flags)
  Siempre escribe ovito_out/counts_<tasa>.csv (clase, color, N).

Entrada : ovito_in/<tasa>/glass300.dump      (id type x y z; 1=Zr 2=Cu 3=Ag)
          ovito_in/<tasa>/glass300.boo.dump  (columnas c_cQ12[1..4] = Q4 Q6 W4hat W6hat, nnn=12)
Salida  : ovito_out/<tasa>_snapshot.png  (o <tasa>_bare.png con --bare) y ovito_out/counts_<tasa>.csv

Clases (mutuamente excluyentes, gana la primera que cumple):
    1 ico-BOO      Q6 >= 0.6 y W6hat <= -0.16                (BOO, nnn=12)
    2 ico-Voronoi  indice exacto <0,0,12,0>
    3 Z16-Voronoi  indice exacto <0,0,12,4>
    4 Z15-Voronoi  indice exacto <0,0,12,3>
    5 Z14-Voronoi  indice exacto <0,0,12,2>
    6 resto
Indice exacto <0,0,12,k> = 12 caras pentagonales, k hexagonales y NINGUNA otra (n3=n4=0, n>=7=0).
Nunca se usa CN==14/15/16 para clasificar; el CN solo se imprime como chequeo del paso Voronoi.
Voronoi: VoronoiAnalysisModifier(compute_indices=True), SIN otros argumentos (default puro).
"""
import sys
from pathlib import Path
import numpy as np
from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier, SliceModifier, ExpressionSelectionModifier, DeleteSelectedModifier
from ovito.data import ParticleType
from ovito.vis import Viewport, TachyonRenderer, OpenGLRenderer, TextLabelOverlay
from ovito.qt_compat import QtCore

IN_DIR = Path("ovito_in")
OUT_DIR = Path("ovito_out")
ALL_RATES = ["2e10", "1e11", "1e12"]
SPECIES = {1: "Zr", 2: "Cu", 3: "Ag"}
Q6_MIN, W6_MAX = 0.6, -0.16
IMG_SIZE = (1600, 1600)
USE_ICO_VOR = True          # --no-ico-voronoi lo pone en False (y quita la clase 2 de CLASSES)
TAGV = ""                   # "noicovor_" en la variante: se antepone a los nombres de salida
CENTER_T = None             # --center <especie>: id de tipo (1=Zr 2=Cu 3=Ag) cuyos centros se dibujan; None = todas
HIDE_RESTO = False          # --resto-transp 1: elimina los atomos "resto" del render (solo visual)
# ico-BOO a 300 K en la Etapa 3 (tabla_clusters_51.csv, Ag+Cu+Zr): chequeo cruzado informativo
EXPECT_ICO_BOO = {"2e10": 14 + 755 + 1, "1e11": 8 + 567 + 0, "1e12": 9 + 428 + 1}

# id, nombre, color RGB, radio (A), transparencia
CLASSES = [
    (1, "ico-BOO",     (0.80, 0.00, 0.00), 1.60, 0.0),
    (2, "ico-Voronoi", (0.95, 0.40, 0.00), 1.50, 0.0),
    (3, "Z16-Voronoi", (0.60, 0.20, 0.80), 1.35, 0.0),
    (4, "Z15-Voronoi", (0.10, 0.35, 0.95), 1.35, 0.0),
    (5, "Z14-Voronoi", (0.05, 0.70, 0.55), 1.35, 0.0),
    (6, "resto",       (0.62, 0.62, 0.62), 0.45, 0.5),
]


# ----------------------------------------------------------------------------- BOO
def read_boo(path):
    """Parsea a mano el dump BOO. Devuelve ids, type, q6, w6hat, longitudes de caja.
    Las columnas se buscan POR NOMBRE en la linea 'ITEM: ATOMS' (c_cQ12[2]=Q6, c_cQ12[4]=W6hat)."""
    with open(path) as f:
        assert f.readline().startswith("ITEM: TIMESTEP")
        f.readline()
        assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
        n = int(f.readline())
        assert f.readline().startswith("ITEM: BOX BOUNDS pp pp pp")
        L = []
        for _ in range(3):
            lo, hi = f.readline().split()[:2]
            L.append(float(hi) - float(lo))
        header = f.readline().split()
        assert header[:2] == ["ITEM:", "ATOMS"]
        cols = header[2:]
        arr = np.loadtxt(f)
    assert arr.shape == (n, len(cols)), f"{path}: {arr.shape} vs {n} x {len(cols)}"
    ix = {c: cols.index(c) for c in ("id", "type", "c_cQ12[2]", "c_cQ12[4]")}   # KeyError/ValueError si faltan
    q6, w6 = arr[:, ix["c_cQ12[2]"]], arr[:, ix["c_cQ12[4]"]]
    assert 0.0 <= q6.min() and q6.max() <= 1.0, f"Q6 fuera de [0,1]: {q6.min()}..{q6.max()}"
    assert -1.0 <= w6.min() and w6.max() <= 1.0, f"W6hat fuera de rango: {w6.min()}..{w6.max()}"
    return (arr[:, ix["id"]].astype(np.int64), arr[:, ix["type"]].astype(np.int64), q6, w6, np.array(L))


def align_boo(ids_pos, types_pos, cell_L, boo):
    """Reordena Q6/W6hat al orden de particulas del pipeline (por id). Asserts explicitos."""
    ids_b, types_b, q6, w6, L_b = boo
    n = len(ids_pos)
    assert len(ids_b) == n, f"N distinto: posiciones {n} vs BOO {len(ids_b)}"
    assert len(np.unique(ids_pos)) == n and len(np.unique(ids_b)) == n, "ids repetidos"
    assert set(ids_pos.tolist()) == set(ids_b.tolist()), "los conjuntos de id difieren"
    assert np.allclose(cell_L, L_b, atol=1e-6), f"caja distinta: {cell_L} vs {L_b}"
    order = np.argsort(ids_b)
    pos_in_sorted = np.searchsorted(ids_b[order], ids_pos)
    idx = order[pos_in_sorted]                      # idx[k] = fila del BOO con el id de la particula k
    assert (ids_b[idx] == ids_pos).all()
    assert (types_b[idx] == types_pos).all(), "type distinto entre dump de posiciones y dump BOO"
    return q6[idx], w6[idx]


# ----------------------------------------------------------------------------- Voronoi
def edge_column_offset(vi):
    """OVITO 'Voronoi Index': columna j = numero de caras con (j + o) aristas. Se determina o por la
    identidad de Euler de celdas simples (3 caras por vertice): sum_k k*n_k = 6*(F-2), F = sum n_k."""
    F = vi.sum(axis=1)
    best = None
    for o in (0, 1, 2):
        edges = np.arange(vi.shape[1]) + o
        frac = np.mean((vi * edges).sum(axis=1) == 6 * (F - 2))
        if best is None or frac > best[1]:
            best = (o, frac)
    o, frac = best
    assert frac >= 0.99, f"ninguna convencion de columnas cumple Euler (mejor o={o}, frac={frac:.4f})"
    return o, frac


def exact_masks(vi):
    """Mascaras de indice EXACTO <0,0,12,k> sobre el vector completo de caras."""
    o, frac = edge_column_offset(vi)
    c5, c6 = 5 - o, 6 - o                            # columnas de pentagonos y hexagonos
    F = vi.sum(axis=1)
    def exact(k):
        return (vi[:, c5] == 12) & (vi[:, c6] == k) & (F == 12 + k)   # F == 12+k => n3=n4=0 y n>=7=0
    return {"ico": exact(0), "z16": exact(4), "z15": exact(3), "z14": exact(2)}, o, frac


def cn_check(cn, species):
    """Imprime CN global y por especie; devuelve True si coincide con el patron conocido."""
    print(f"  CN global: media={cn.mean():.3f}")
    hist = np.bincount(cn, minlength=20)
    print("  histograma CN (global): " + " ".join(f"{k}:{hist[k]}" for k in range(8, 20) if hist[k]))
    ok = 13.5 <= cn.mean() <= 14.5
    modes = {}
    for t, name in SPECIES.items():
        h = np.bincount(cn[species == t], minlength=20)
        modes[name] = int(np.argmax(h))
        print(f"  CN {name}: media={cn[species == t].mean():.2f} moda={modes[name]}")
    ok &= (modes["Cu"] == 12 and modes["Ag"] == 14 and modes["Zr"] == 15)
    bimodal = hist[12] > hist[11] and hist[12] > hist[13] and hist[15] > hist[14] and hist[15] > hist[16]
    print(f"  bimodal global en 12 y 15 (maximos locales): {bimodal}")
    print(f"  chequeo CN (media 13.5-14.5, moda Cu=12 / Ag=14 / Zr=15): {'OK' if ok else 'FALLA'}")
    return ok


# ----------------------------------------------------------------------------- clasificacion
def classify(boo_ico, masks):
    """Prioridad ico-BOO > ico-Vor > Z16 > Z15 > Z14 > resto. Devuelve clase (1..6) por particula.
    Con USE_ICO_VOR=False se salta ico-Vor (esos atomos quedan en resto salvo que ya sean ico-BOO)."""
    cls = np.full(len(boo_ico), 6, dtype=np.int32)
    taken = np.zeros(len(boo_ico), dtype=bool)
    order = [(1, boo_ico)] + ([(2, masks["ico"])] if USE_ICO_VOR else []) + \
            [(3, masks["z16"]), (4, masks["z15"]), (5, masks["z14"])]
    for cid, m in order:
        sel = m & ~taken
        cls[sel] = cid
        taken |= sel
    return cls


def report(rate, cls, species, boo_ico, masks):
    names = {c[0]: c[1] for c in CLASSES}
    print(f"  --- clases excluyentes (las del render) [{rate}] ---")
    print(f"  {'clase':<12}{'N':>7}" + "".join(f"{n:>7}" for n in SPECIES.values()))
    for cid in names:
        m = cls == cid
        print(f"  {names[cid]:<12}{m.sum():>7}" + "".join(f"{(m & (species == t)).sum():>7}" for t in SPECIES))
    print(f"  suma = {np.bincount(cls, minlength=7)[1:].sum()}")
    print("  --- mascaras crudas SIN prioridad (categorias distintas, no se suman) ---")
    for label, m in (("ico-BOO", boo_ico), ("ico-Voronoi", masks["ico"]), ("Z16-Voronoi", masks["z16"]),
                     ("Z15-Voronoi", masks["z15"]), ("Z14-Voronoi", masks["z14"])):
        print(f"  {label:<12}{m.sum():>7}" + "".join(f"{(m & (species == t)).sum():>7}" for t in SPECIES))
    exp = EXPECT_ICO_BOO[rate]
    got = int(boo_ico.sum())
    print(f"  ico-BOO vs Etapa 3 (tabla_clusters_51, 300 K): {got} vs {exp} -> {'OK' if got == exp else 'DIFERENTE'}")


# ----------------------------------------------------------------------------- render
def make_assign_modifier(cls):
    def assign(frame, data):
        del data.particles_["Particle Type"]      # si no, se reusan los tipos 1,2,3 del archivo (ids duplicados)
        types = data.particles_.create_property("Particle Type", data=cls)
        for cid, name, color, radius, _ in CLASSES:
            types.types.append(ParticleType(id=cid, name=name, color=color, radius=radius))
        assert [t.id for t in types.types] == [c[0] for c in CLASSES], "tipos duplicados o faltantes"
        transp = np.zeros(len(cls), dtype=np.float32)
        for cid, _, _, _, tr in CLASSES:
            transp[cls == cid] = tr
        data.particles_.create_property("Transparency", data=transp)
    return assign


def render(pipeline, rate, cls, out_png, opengl=False, bare=False, slab=None):
    counts = np.bincount(cls, minlength=7)
    pipeline.add_to_scene()
    try:
        pipeline.source.data.cell_.vis.line_width = 0.6
        if slab:
            vp = Viewport(type=Viewport.Type.Top)                    # ortografica, mirando a -z
            vp.zoom_all(size=IMG_SIZE)
        else:
            vp = Viewport(type=Viewport.Type.Perspective, camera_dir=(-1.0, -1.25, -0.8))
            vp.zoom_all(size=IMG_SIZE)
            vp.fov = 2 * np.arctan(1.22 * np.tan(vp.fov / 2))        # zoom-out ~20%: deja margen para la leyenda
        rate_txt = {"2e10": "2×10¹⁰", "1e11": "1×10¹¹", "1e12": "1×10¹²"}[rate]
        AF = QtCore.Qt.AlignmentFlag

        def label(text, color, align, x, y, size):
            vp.overlays.append(TextLabelOverlay(text=text, alignment=align, offset_x=x, offset_y=y,
                                                font_size=size, text_color=color,
                                                outline_enabled=True, outline_color=(1, 1, 1)))
        if not bare:                                                  # --bare: sin texto (para el compuesto)
            label(f"{rate_txt} K/s, 300 K" + (f" (slab {slab:g} Å)" if slab else ""), (0, 0, 0), AF.AlignLeft | AF.AlignTop, 0.02, -0.02, 0.03)
            for i, (cid, name, color, _, _) in enumerate(CLASSES):   # leyenda: 2 filas x 3 columnas, abajo
                label(f"● {name}: {counts[cid]}", color, AF.AlignLeft | AF.AlignBottom,
                      0.02 + 0.325 * (i % 3), 0.02 + 0.04 * (1 - i // 3), 0.02)
        OUT_DIR.mkdir(exist_ok=True)
        renderer = OpenGLRenderer() if opengl else TachyonRenderer(direct_light_intensity=0.7,
                                                                   ambient_occlusion_brightness=0.7)
        vp.render_image(size=IMG_SIZE, filename=str(out_png), background=(1, 1, 1), renderer=renderer)
    finally:
        pipeline.remove_from_scene()
    print(f"  -> {out_png}")


# ----------------------------------------------------------------------------- main
def write_counts(rate, cls, tag="", n_total=None):
    """counts_<tasa>.csv: fuente unica de colores/nombres/conteos para compose_snapshots.py."""
    OUT_DIR.mkdir(exist_ok=True)
    counts = np.bincount(cls, minlength=7)
    with open(OUT_DIR / f"counts_{tag}{TAGV}{rate}.csv", "w") as f:
        f.write("rate,class_id,name,r,g,b,N,n_total\n")
        for cid, name, color, _, _ in CLASSES:
            f.write(f"{rate},{cid},{name},{color[0]},{color[1]},{color[2]},{counts[cid]},{n_total or len(cls)}\n")


def process(rate, do_render, force, opengl, bare, slab):
    print(f"=== {rate} ===")
    pos_path, boo_path = IN_DIR / rate / "glass300.dump", IN_DIR / rate / "glass300.boo.dump"
    pipeline = import_file(str(pos_path))
    pipeline.modifiers.append(VoronoiAnalysisModifier(compute_indices=True))     # default puro
    data = pipeline.compute()
    ids = np.asarray(data.particles["Particle Identifier"][...], dtype=np.int64)
    species = np.asarray(data.particles["Particle Type"][...], dtype=np.int64)
    cell_L = np.diag(np.asarray(data.cell[...])[:, :3])
    vi = np.asarray(data.particles["Voronoi Index"][...], dtype=np.int64)
    cn = np.asarray(data.particles["Coordination"][...], dtype=np.int64)
    assert set(np.unique(species)) == {1, 2, 3} and len(ids) == 65536
    assert (vi.sum(axis=1) == cn).all(), "Voronoi Index no suma la Coordination"

    q6, w6 = align_boo(ids, species, cell_L, read_boo(boo_path))
    boo_ico = (q6 >= Q6_MIN) & (w6 <= W6_MAX)

    masks, o, frac = exact_masks(vi)
    print(f"  Voronoi Index: columna j = caras con (j+{o}) aristas (Euler ok en {100 * frac:.2f}% de los atomos)")
    cn_ok = cn_check(cn, species)
    cls = classify(boo_ico, masks)
    report(rate, cls, species, boo_ico, masks)
    write_counts(rate, cls)
    if CENTER_T:                                                     # solo visual: otras especies -> "resto" (oculto/gris)
        of_sp = species == CENTER_T
        cls = np.where(of_sp, cls, 6).astype(np.int32)
        write_counts(rate, cls[of_sp], n_total=int(of_sp.sum()))     # N por clase DENTRO de la especie, % sobre esa especie
        print(f"  --center {SPECIES[CENTER_T]}: {int(of_sp.sum())} atomos de {SPECIES[CENTER_T]}; clases (N): "
              + ", ".join(f"{c[1]}={int((cls == c[0]).sum())}" for c in CLASSES if c[0] != 6))

    if not cn_ok and not force:
        print("  !! chequeo de CN fallo -> no se renderiza (usar --force para seguir)")
        return False
    if do_render:
        pipeline.modifiers.append(make_assign_modifier(cls))
        tp = pipeline.compute().particles["Particle Type"]
        for t, (cid, name, color, radius, _) in zip(tp.types, CLASSES):
            assert t.id == cid and t.name == name and np.allclose(t.color, color) and abs(t.radius - radius) < 1e-9
        assert (np.asarray(tp[...]) == cls).all()
        cls_show = cls
        if slab:                                                     # corte SOLO visual, tras clasificar todo el sistema
            zc = float(np.asarray(data.cell[...])[2, 3]) + cell_L[2] / 2
            pipeline.modifiers.append(SliceModifier(normal=(0, 0, 1), distance=zc, slab_width=slab))
            cls_show = np.asarray(pipeline.compute().particles["Particle Type"][...], dtype=np.int32)
            print(f"  slab {slab:g} A centrado en z={zc:.2f}: {len(cls_show)} atomos visibles "
                  f"({100 * len(cls_show) / len(cls):.1f}% del sistema)")
            write_counts(rate, cls_show, tag="slab_")
        if HIDE_RESTO:                                               # solo visual, DESPUES de contar/cortar
            pipeline.modifiers.append(ExpressionSelectionModifier(expression="ParticleType==6"))
            pipeline.modifiers.append(DeleteSelectedModifier())
        name = ("slab_" if slab else "") + TAGV + ("bare" if bare else "snapshot")
        render(pipeline, rate, cls_show, OUT_DIR / f"{rate}_{name}.png", opengl, bare, slab)
    return True


if __name__ == "__main__":
    args = sys.argv[1:]
    slab = None
    if "--slab" in args:
        i = args.index("--slab")
        slab = float(args[i + 1])
        del args[i:i + 2]
    if "--no-ico-voronoi" in args:
        USE_ICO_VOR, TAGV = False, "noicovor_"
        CLASSES = [c for c in CLASSES if c[0] != 2]
    if "--resto-transp" in args or "--radius-scale" in args:
        rt, rs = 0.5, 1.0                                            # defaults = version original
        if "--resto-transp" in args:
            i = args.index("--resto-transp"); rt = float(args[i + 1]); del args[i:i + 2]
        if "--radius-scale" in args:
            i = args.index("--radius-scale"); rs = float(args[i + 1]); del args[i:i + 2]
        assert 0.0 <= rt <= 1.0 and rs > 0
        HIDE_RESTO = rt >= 1.0
        CLASSES = [(cid, n, col, r if cid == 6 else r * rs, rt if cid == 6 else tr) for cid, n, col, r, tr in CLASSES]
        TAGV += f"rt{rt:g}_rs{rs:g}_"
    if "--center" in args:
        i = args.index("--center"); sp_name = args[i + 1]; del args[i:i + 2]
        assert sp_name in SPECIES.values(), f"especies validas: {list(SPECIES.values())}"
        CENTER_T = {v: k for k, v in SPECIES.items()}[sp_name]
        TAGV += f"{sp_name}_"
    do_render, force, opengl, bare = ("--no-render" not in args, "--force" in args, "--opengl" in args, "--bare" in args)
    rates = [a for a in args if not a.startswith("--")] or ALL_RATES
    assert all(r in ALL_RATES for r in rates), f"tasas validas: {ALL_RATES}"
    results = [process(r, do_render, force, opengl, bare, slab) for r in rates]
    sys.exit(0 if all(results) else 1)
