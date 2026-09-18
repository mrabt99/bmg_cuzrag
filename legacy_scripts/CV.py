#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 09:02:58 2023

@author: sebas
"""

import numpy as np
import lammps_logfile
from lammps_logfile import running_mean
import matplotlib.pyplot as plt

path = '/home/sebas/Documentos/BMG/corridas/4754755/cooling/10/'

log = lammps_logfile.File(path + "log.lammps")

key = log.get_keywords()
print("Log keywords:", key)

nparticles = 8192
step = log.get(key[10], run_num=1)
temp = log.get(key[11], run_num=1)
tot_eng = log.get(key[12], run_num=1)
pot_eng = log.get(key[2], run_num=1)
entalpy = log.get(key[0], run_num=1)


n_mean = 100  # default 50
temp_ave = running_mean(temp,n_mean)
ene_ave = running_mean(tot_eng, n_mean)
pot_ave = running_mean(pot_eng, n_mean)
ental_ave = running_mean(entalpy, n_mean)


x = []
y = []

for i, values in enumerate(temp_ave):
    if np.isnan(values) == False:
        y.append(ene_ave[i])

# npuntos = 200
npuntos = 50

t_ini = 2000-(n_mean/2)
t_fin = 300+(n_mean/2)

x = np.linspace(t_ini,t_fin,npuntos)


paso = int(len(y)/npuntos)

ene = []
for i in range(0,len(y)-paso,paso):
    ene.append(y[i])
ene = np.array(ene)

fig, ax = plt.subplots()
ax.scatter(x,ene, s=10)
#ax.scatter(x,ene,'o')



dt = x[0] - x[1]

# Cp a partir de la energia
cp = [(ene[1]-ene[0])/dt]
for i in range(len(x)-1):
    cp.append((ene[i-1]-ene[i])/dt)
    
plt.plot(x,cp,'*')
plt.ylim(2,3.5)
plt.show()


x_ave = running_mean(x,5)
y_ave = running_mean(cp,5)
plt.plot(x_ave,y_ave,'*')
plt.ylim(2,3.5)
plt.show()

####################################################
y = []

for i, values in enumerate(temp_ave):
    if np.isnan(values) == False:
        y.append(pot_ave[i])

pot = []
for i in range(0,len(y)-paso,paso):
    pot.append(y[i])
pot = np.array(pot)

#fig, ax = plt.subplots()
#ax.scatter(x,pot, s=10)

dt = x[0] - x[1]

# Cp a partir de la energia
cp2 = [(pot[1]-pot[0])/dt]
for i in range(len(x)-1):
    cp2.append((pot[i-1]-pot[i])/dt)
    
plt.plot(x,cp2,'r*')
plt.ylim(1,2.5)
plt.show()


x_ave = running_mean(x,5)
y_ave = running_mean(cp2,5)
plt.plot(x_ave,y_ave,'r*')
plt.ylim(1,2.5)
plt.show()


####################################################


# Cp a partir de la entalpia

y_ental = []

for i, values in enumerate(temp_ave):
    if np.isnan(values) == False:
        y_ental.append(ental_ave[i])
         
cp_ental = [(y_ental[1]-y_ental[0])/dt]
for i in range(len(x)-1):
    cp_ental.append((y_ental[i-1]-y_ental[i])/dt)
    
plt.plot(x,cp_ental,'*')
plt.ylim(-2,2)
plt.show()


x_ave = running_mean(x,20)
y_ental_ave = running_mean(cp_ental,20)
plt.plot(x_ave,y_ental_ave,'*')
plt.ylim(0,1)

