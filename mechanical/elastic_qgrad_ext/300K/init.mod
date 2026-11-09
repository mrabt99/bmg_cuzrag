# ======================================================================
#  ELASTICAS A 300 K  (adaptado de examples/ELASTIC_T/DEFORMATION, LAMMPS stable_23Jun2022)
#  Parametros: se cambian con  -var nombre valor  en la linea de comandos.
# ======================================================================
variable up      index 0.005       # deformacion (ingenieria) de cada perturbacion: +-up
variable temp    index 300.0       # K
variable tdamp   index 0.1         # ps, termostato Nose-Hoover (NVT)
variable pdamp   index 1.0         # ps, barostato (solo el NPT de preparacion)
variable nnpt    index 50000       # pasos del NPT que fija la caja PROMEDIO (50 ps)
variable nequil  index 20000       # pasos de equilibrado tras cada perturbacion (20 ps)
variable nrun    index 100000      # pasos de promedio de las tensiones (100 ps)
variable nser    index 1000        # ventana (pasos) de los promedios parciales, para estimar el error (1 ps)
variable rst     index /scratch/negishi/aaabt/bmg_cuzrag/qgrad_ext/restart.glass300K_N65536_gradext
variable pot     index /scratch/negishi/aaabt/bmg_cuzrag/rate_2e10/ZrCuAg.lammps.eam

# unidades metal: la presion sale en bar -> constantes elasticas en GPa (1 bar = 1e-4 GPa)
variable cfac    equal 1.0e-4
variable cunits  string GPa
variable nevery  equal 10                     # intervalo de muestreo de fix ave/time
# fix ave/time cierra sus ventanas en pasos multiplos de Nfreq (absolutos): por eso los promedios de la caja y de la presion usan
# ventanas de nser pasos con 'ave running' (promedio de todas las ventanas completas desde que se define el fix), no una sola ventana de nrun pasos.
variable nreps   equal ${nser}/${nevery}
