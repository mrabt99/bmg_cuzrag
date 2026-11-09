# Potencial y salida: IGUAL que el temple y la traccion (pair_coeff copiado tal cual; tipos 1=Zr 2=Cu 3=Ag)
pair_style      eam/alloy
pair_coeff      * * ${pot} Zr Cu Ag
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes
timestep        0.001
thermo          ${nser}
thermo_style    custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol density
thermo_modify   norm no
