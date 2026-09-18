#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:02:41 2023

@author: sebas
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Times'
rcParams['font.size'] = 10
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['xtick.direction'] = 'in'
rcParams['xtick.major.size'] = 3
rcParams['xtick.major.width'] = 0.5
rcParams['xtick.minor.size'] = 1.5
rcParams['xtick.minor.width'] = 0.5
rcParams['xtick.top'] = True
rcParams['ytick.direction'] = 'in'
rcParams['ytick.major.size'] = 3
rcParams['ytick.major.width'] = 0.5
rcParams['ytick.minor.size'] = 1.5
rcParams['ytick.minor.width'] = 0.5
rcParams['ytick.right'] = True
rcParams['lines.linewidth'] = 1.
rcParams['axes.linewidth'] = 0.5

path = '/home/sebas/Documentos/BMG/corridas/4754755/cooling/10/boo/'
z = np.loadtxt(path+'dump.110000000.cc')
zr = 1
cu = 2
ag = 3
ztypes = z[:,1]
q4z = []
q6z = []
q4c = []
q6c = []
q4a = []
q6a = []

for i, ztype  in enumerate(ztypes):
    if ztype == zr:
        q4z.append(z[i,2])
        q6z.append(z[i,3])
    if ztype == cu:
        q4c.append(z[i,2])
        q6c.append(z[i,3])
    if ztype == ag:
        q4a.append(z[i,2])
        q6a.append(z[i,3])   

plt.plot(q4z,q6z,'b*')
plt.plot(q4c,q6c,'g*')
plt.plot(q4a,q6a,'r*')
plt.show()

num_bins = 100
n, bins, patches = plt.hist(q6z, num_bins, facecolor ='blue')
n1, bins1, patches1 = plt.hist(q6c, num_bins, facecolor ='green')
n2, bins2, patches2 = plt.hist(q6a, num_bins, facecolor ='red')
plt.plot(bins,bins1)
plt.show()


