# ======================================================================
#  ELASTICAS A 0 K  (adaptado de examples/ELASTIC, LAMMPS stable_23Jun2022: minimizacion tras cada perturbacion)
#  Parametros: se cambian con  -var nombre valor.
# ======================================================================
variable up      index 1.0e-4      # deformacion (ingenieria) de cada perturbacion: +-up  (el ejemplo usa 1e-6 para cristales; un vidrio necesita mas)
variable ftol    index 1.0e-8      # tolerancia de fuerza de la minimizacion (eV/A)
variable maxiter index 100000      # iteraciones maximas por minimizacion
variable maxeval index 200000      # evaluaciones de fuerza maximas por minimizacion
variable etol    equal 0.0
variable dmax    equal 1.0e-2
variable rst     index /scratch/negishi/aaabt/bmg_cuzrag/structural_analysis/qgrad/relax300K/restart.glass300K_N65536_qgrad
variable pot     index /scratch/negishi/aaabt/bmg_cuzrag/rate_2e10/ZrCuAg.lammps.eam
# unidades metal: presion en bar -> GPa
variable cfac    equal 1.0e-4
variable cunits  string GPa
