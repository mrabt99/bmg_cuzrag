# Potencial IGUAL que el temple y la traccion (pair_coeff copiado tal cual; tipos 1=Zr 2=Cu 3=Ag) + estilo de minimizacion del ejemplo
pair_style      eam/alloy
pair_coeff      * * ${pot} Zr Cu Ag
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes
min_style       cg
min_modify      dmax ${dmax} line quadratic
thermo          1000
thermo_style    custom step pe press pxx pyy pzz pxy pxz pyz lx ly lz vol density fmax fnorm
thermo_modify   norm no
