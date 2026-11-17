# Prueba: `stress/atom` contra `centroid/stress/atom` en `heat/flux` (16/11/2026)

**Motivo.** Se propuso cambiar `compute stress/atom` por `compute centroid/stress/atom`. El argumento es que con potenciales de muchos cuerpos el reparto de a pares del virial daría un flujo de calor espurio.

**Predicción.** En EAM la fuerza entre i y j es central y de a pares, F_ij = −[φ′ + F′_i ρ′_j + F′_j ρ′_i] r̂_ij. El reparto del virial a medias coincide con el del centroide, así que J tiene que salir igual. El centroide cambia el resultado solo con términos angulares (Tersoff, SW, MEAM, ángulos o diedros).

**Qué se corrió.** Job 46000627, `-q standby`, 48 tareas, 34 s. `in.centroid` lee `../cal_300K_2e10/cal_npt.restart`, reescalado a `vmean`. Pone velocidades nuevas a 300 K con semilla 1001, termaliza 2 ps en NVT y corre 5 ps en NVE. En la NVE calcula a la vez los dos `heat/flux` sobre la misma trayectoria, uno con cada tensor, cada 5 pasos. También guarda el máximo, sobre todos los átomos, de |S_st − S_cst| y de la asimetría xy − yx del centroide.

**Resultado** (`compara.py` → `compara.txt`):
- J total y J convectivo: max|ΔJ| = **0,0 exacto** en las 1001 muestras, con rms(J) = 563 eV·Å/ps. ⟨J²⟩/3 da 3,166468e5 en los dos casos.
- Tensor por átomo: max|S_st − S_cst| = 0,0 en las 6 componentes. La asimetría del centroide también es 0,0, y la escala del tensor es 1,2e7 bar·Å³.
- Control: la corrida está sana. T ≈ 300 K y E_tot = −314473,56 eV constante durante toda la NVE.

**Lectura.** La igualdad es bit a bit, porque LAMMPS marca `eam/alloy` como `CENTROID_SAME` y usa el mismo virial de pares en los dos computes. Cambiar el compute no cambia κ. La corrección física de la fórmula de pares en EAM la dan el argumento de las fuerzas centrales y el acuerdo de GK (0,78) con NEMD (0,73, `../nemd/`), que no usa el virial por átomo.
