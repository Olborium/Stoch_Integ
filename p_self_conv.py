###########################
# Strong convergence test #
###########################
#
# The script reads M groups of trajectories computed using M reference Brownian processes with time step dt.
# Each group contains N trajectories computed with N different time steps DT >= dt.
# All trajectories are computed using the same scheme.
# For each group the script finds the absolute difference dist(DT) between the end points of the trajectories with time steps dt and DT.
# The averaging is performed over M groups producing dist_av(DT) and statistical uncertainty var(DT).
# The data points dist_av(DT) are fitted by a line (slope a, intercept b) in the log-log scale.
# The measured strong order of convergence equals a.
# The script prints the pair (a, b) and makes a plot showing dist_av(DT) and the linear fit.
#
# The stochastic pendulum is used in this example script.

import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})

TEMP = float(sys.argv[1])
ETA = float(sys.argv[2])
SIMTIME = float(sys.argv[3])
BETA = int(sys.argv[4])
GAMMA = int(sys.argv[5])
SAMPLES = int(sys.argv[6])
INTEG = sys.argv[7]
M = int(sys.argv[8])
N = int(sys.argv[9])
Ns = np.zeros(N)
for i in range(N):
    Ns[i] = int(sys.argv[i+10])

path = 'Out/Pendulum/SelfConv/'+INTEG
files = os.listdir(path)
Qfin = np.zeros((N,M))
Pfin = np.zeros((N,M))
str = f'{Ns[0]:.0f}'+'_'+f'{BETA:.0f}'+'_'+f'{GAMMA:.0f}'
for filename in files:
    for j in range(M):
        for i in range(N):
            str = f'{Ns[i]:.0f}'+'_'+f'{BETA:.0f}'+'_'+f'{GAMMA:.0f}'+'_'+f'{j:.0f}'                     
            if re.match('pendulum_qt_'+str,filename):
                filepath = os.path.join(path, filename)
                with open(filepath, 'r') as file:
                    line = file.readline().strip()
                    elements = line.split()
                    last_element = elements[1]
                    if np.isnan(float(last_element)):
                        continue                
                Qfin[i,j] = np.loadtxt(filepath)[-1]    
for filename in files:
    for j in range(M):
        for i in range(N):
            str = f'{Ns[i]:.0f}'+'_'+f'{BETA:.0f}'+'_'+f'{GAMMA:.0f}'+'_'+f'{j:.0f}'                
            if re.match('pendulum_pt_'+str,filename):
                filepath = os.path.join(path, filename)
                with open(filepath, 'r') as file:
                    line = file.readline().strip()
                    elements = line.split()
                    last_element = elements[1]
                    if np.isnan(float(last_element)):
                        continue                
                Pfin[i,j] = np.loadtxt(filepath)[-1]            

dist = np.zeros((N-1,M))
for i in range(N-1):
    for j in range(M):
        dist[i,j] = np.sqrt((Qfin[i+1,j]-Qfin[0,j])**2+(Pfin[i+1,j]-Pfin[0,j])**2)

dist_av = np.zeros(N-1)
var = np.zeros(N-1)
for i in range(N-1):
    dist_av[i] = np.sum(dist[i,:])/M
    var[i] = np.sqrt(np.sum((dist[i,:]-dist_av[i])**2))/(M-1)

def f(h,a,b):
    return a*h + b

DTs = np.zeros(N)
for i in range(N):
    DTs[i] = SIMTIME/Ns[i]

iMax = np.where(dist_av<1)[0][-1]+1
params, cov = curve_fit(f, np.log10(DTs[1:iMax+1]), np.log10(dist_av[:iMax]))
(a, b) = (params[0], params[1])
print(a,b)

t = np.linspace(DTs[1],DTs[-1],100)
str = f'{BETA}'+f'{GAMMA}'+'_M='+f'{M:.0f}'+'_dt='+f'{DTs[0]:.5f}'+'_eta='+f'{ETA:.2f}'+'_int_'+INTEG+'x'
if b>0:
    str1 = r'$y=$'+f'{a:.2f}'+r'$x+$'+f'{np.abs(b):.2f}'
else:
    str1 = r'$y=$'+f'{a:.2f}'+r'$x-$'+f'{np.abs(b):.2f}'
fig, ax = plt.subplots(figsize=(6,4))
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=14)
ax.set_xlabel(r'$h$',fontsize=18)
ax.set_ylabel(r'$\langle | \Delta z_T | \rangle$',fontsize=18)
ax.set_xscale('log',base=10)
ax.set_yscale('log',base=10)
plt.ylim((0.5*dist_av[0],50))
ax.plot(t, t**a*10**b,'blue',linestyle='dashed')
ax.errorbar(DTs[1:], dist_av, var, color='black', fmt='o', ms=5, markerfacecolor='black', elinewidth=1, capsize=2)
plt.text(DTs[1],dist_av[iMax-1],str1,fontsize=14,color='blue')
plt.savefig('Plots/pend_self_'+str+'.pdf', bbox_inches='tight')