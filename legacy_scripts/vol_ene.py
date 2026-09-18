# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import lammps_logfile
from lammps_logfile import running_mean
import matplotlib.pyplot as plt

path = '/media/sebas/Mis_archivos/BMG/corridas/4754755/cooling/10/'

log = lammps_logfile.File(path + "log.lammps")

key = log.get_keywords()
print("Log keywords:", key)

nparticles = 8192
step = log.get(key[10], run_num=1)
temp = log.get(key[11], run_num=1)
pot_eng = log.get(key[2], run_num=1)
vol = log.get(key[13], run_num=1)
ental = log.get(key[0], run_num=1)

n_mean = 50
temp_ave = running_mean(temp,n_mean)
vol_ave = running_mean(vol, n_mean)
ene_ave = running_mean(pot_eng, n_mean)
h_ave = running_mean(ental, n_mean)



x = []
y = []
e = []
h = []
len_ave = 60
len_ave2 = int(len_ave/2)
for i in range(len_ave2,len(temp_ave)-len_ave2):
    if int(temp_ave[i])%len_ave == 0:
        x.append(temp_ave[i])
        y.append(vol_ave[i])
        e.append(ene_ave[i]/nparticles)
        h.append(h_ave[i])

plt.plot(x,np.array(y)/1000,'*')
plt.show()

x = np.array(x)
y = np.array(y)
e = np.array(e)
h = np.array(h)
x_high = x[:len(x)-9]
y_high = y[:len(y)-9]/1000
x_low = x[len(x)-7:]
y_low = y[len(x)-7:]/1000

A = np.vstack([x_high, np.ones(len(x_high))]).T
m_high, c_high = np.linalg.lstsq(A,y_high,rcond=None)[0]
B = np.vstack([x_low, np.ones(len(x_low))]).T
m_low, c_low = np.linalg.lstsq(B,y_low,rcond=None)[0]

#plt.plot(x_high,y_high,'*')
#plt.plot(x_low,y_low,'*')
plt.plot(x,y/1000,'*')
plt.plot(x_high, m_high*x_high + c_high, 'b')
plt.plot(x_low, m_low*x_low + c_low)    
plt.show()

xp = np.linspace(638.5,638.6,1000)
plt.plot(xp, m_high*xp + c_high)
plt.plot(xp, m_low*xp + c_low)
plt.show()

plt.plot(x,e,'*')
plt.show()





