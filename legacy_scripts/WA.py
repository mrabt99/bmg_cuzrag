#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 10:32:33 2023

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
from lammps_logfile import running_mean
from scipy import stats
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

path = '/home/sebas/Documentos/BMG/corridas/4754755/cooling/10/gofr/'
end_dat = 1710
temp_ini = 2000
temp_fin = 300
temp = np.linspace(temp_ini,temp_fin,end_dat)
col = 4
g = []
l = []
m = []
temp = temp.astype(int)
for i in range(end_dat):
   z = np.loadtxt(path + 'gofr4754755.data'+ str(i)+'_gofr.dat')
   g.append(z[:,col].max())
   ind_max  = z[:,col].argmax()
   l.append(z[ind_max,1])
   m.append(z[ind_max:,col].min())
   #ind_min =  z[ind_max:,col].argmin()
   #m.append(z[ind_min,1])



n_mean = 100
g_ave = running_mean(g,n_mean)
m_ave = running_mean(m,n_mean)
l_ave = running_mean(l,n_mean)

x = []
y = []
z = []
w = []
len_ave = 50
len_ave2 = int(len_ave/2)
for i in range(len_ave2,len(temp)-len_ave2):
    if int(temp[i])%len_ave == 0:
        x.append(temp[i])
        y.append(g_ave[i])
        z.append(l_ave[i])
        w.append(m_ave[i])

x = np.array(x)
plt.plot(x,(np.array(w)/np.array(y)),'*')
plt.show()

wa = np.array(w)/np.array(y)
res1 = stats.linregress(x[26:32],wa[26:32])
res2 = stats.linregress(x[:25],wa[:25])

ys = [0.00,0.015,0.030,0.045,0.060,0.075]
fig1, ax1 = plt.subplots()
ax1.plot(x, wa, 'o', mew=2, alpha=0.8)
ax1.set_xlabel(r'$T\mathrm{(K)}$')
ax1.set_ylabel(r'MWA')
ax1.set_yticks(ys)
#ax1.yaxis.set_major_locator(plt.MaxNLocator(6))
ys2 = [0.00,0.01,0.020,0.03,0.04]
ax2 = plt.axes([0,0,1,1])
ip = InsetPosition(ax1, [0.11,0.54,0.4,0.4])
ax2.set_axes_locator(ip)
ax2.plot(x[14:],wa[14:],'o')
ax2.plot(x[20:],res1[0]*x[20:]+res1[1])
ax2.plot(x[14:29],res2[0]*x[14:29]+res2[1])
ax2.set_yticks(ys2)
ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
#plt.savefig('new_mwa.pdf')
