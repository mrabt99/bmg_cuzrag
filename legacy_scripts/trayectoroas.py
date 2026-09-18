#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 13:53:58 2023

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import powerlaw as pl

def distancia(v1,v2):
    dife = v1-v2
    return np.sqrt(np.dot(dife,dife))

def cdf(data):
    sorted_data = np.sort(data)
    p = 1. * np.arange(len(sorted_data)) / float(len(sorted_data) - 1)
    fig = plt.figure()
    fig.suptitle('CDF of data points')
    ax2 = fig.add_subplot(111)
    ax2.loglog(sorted_data, 1-p)
    ax2.set_xlabel('sorted_random_data')
    ax2.set_ylabel('p')

    

path = '/media/sebas/Mis_archivos/BMG/corridas/4754755/cooling/10/'

nzr = 0
ncu = 0
nag = 0

#leemos eo archivo y generamos una lista para cada elemento
coord_zr = []
coord_cu = []
coord_ag = []
with open(path+'Cu475Zr475Ag5.lammpstrj') as file:
    lineas = []
    for line in file:
        lineas.append(line)
        sl = line.split()
        if (len(sl) > 2):
             if sl[1] == '1':
                nzr += 1
                coord_zr.append(sl)
             if sl[1] == '2':
                coord_cu.append(sl)
                ncu += 1
             if sl[1] == '3':
                 coord_ag.append(sl)
                 nag +=1

natoms = int(lineas[3])
nframes = int(len(lineas)/(9+natoms))

nzr =  int(nzr/nframes)
ncu =  int(ncu/nframes)
nag =  int(nag/nframes)
# Una matriz para cada elemento matriz para los elementos
zr = np.zeros((nframes, nzr,5))
cu = np.zeros((nframes, ncu,5))
ag = np.zeros((nframes, nag,5))
   
zr2 = np.zeros((nzr,5))
cu2 = np.zeros((ncu,5))
ag2 = np.zeros((nag,5))

for i in range(nframes):
    for j in range(nzr):
        k = nzr * i + j
        zr2[j] = np.array(coord_zr[k])
    zrsorted = zr2[zr2[:, 0].argsort()]
    for j in range(ncu):
        k = ncu * i + j
        cu2[j] = np.array(coord_cu[k])
    cusorted = cu2[cu2[:, 0].argsort()]
    for j in range(nag):
        k = nag * i + j
        ag2[j] = np.array(coord_ag[k])
    agsorted = ag2[ag2[:, 0].argsort()]
    zr[i] = zrsorted
    cu[i] = cusorted
    ag[i] = agsorted
zr = np.delete(zr,1, axis=2)
cu = np.delete(cu,1, axis=2)
ag = np.delete(ag,1, axis=2)


#Distancias recorridas
zr_dist = np.zeros(nzr)
for i in range(nzr):
    dist = 0
    for j in range(1,len(zr)):
        v1 = np.array([zr[j-1,i,1], zr[j-1,i,2], zr[j-1,i,3]])
        v2 = np.array([zr[j,i,1], zr[j,i,2], zr[j,i,3]])                                     
        dist += distancia(v1,v2)
    zr_dist[i] = dist
    
cu_dist = np.zeros(ncu)
for i in range(ncu):
    dist = 0
    for j in range(1,len(cu)):
        v1 = np.array([cu[j-1,i,1], cu[j-1,i,2], cu[j-1,i,3]])
        v2 = np.array([cu[j,i,1], cu[j,i,2], cu[j,i,3]])                                     
        dist += distancia(v1,v2)
    cu_dist[i] = dist
    
ag_dist = np.zeros(nag)
for i in range(nag):
    dist = 0
    for j in range(1,len(ag)):
        v1 = np.array([ag[j-1,i,1], ag[j-1,i,2], ag[j-1,i,3]])
        v2 = np.array([ag[j,i,1], ag[j,i,2], ag[j,i,3]])                                     
        dist += distancia(v1,v2)
    ag_dist[i] = dist
        

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.scatter(ag[0,:,1],ag[0,:,2],ag[0,:,3])

fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
ax2.plot(ag[:,0,1],ag[:,0,2],ag[:,0,3],'-*')

fig3, ax3 = plt.subplots(subplot_kw={"projection": "3d"})
ax3.scatter(ag[17,:,1],ag[17,:,2],ag[17,:,3])

fig4, ax4 = plt.subplots()
ax4.hist(zr_dist,50,density = True, label = 'Zr')
ax4.hist(cu_dist,50,density = True, label = 'Cu')
ax4.hist(ag_dist,50,density = True, label='Ag')
ax4.legend()

#Cumulative Distribution
xzr, yzr = pl.ccdf(zr_dist)
xcu, ycu = pl.ccdf(cu_dist)
xag, yag = pl.ccdf(ag_dist)


fig5, ax5 = plt.subplots()
ax5.loglog(xzr,yzr, label = 'Zr')
ax5.loglog(xcu,ycu, label = 'Cu')
ax5.loglog(xag,yag,label='Ag')
ax5.legend()



