#%%
import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d 
from scipy.optimize import curve_fit

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})
#%%
""" Parameters of the simulations """

SIZE = 100
N_x = 8192
DX = SIZE/N_x
TimeOfSim = 4000
M = 13
DT = 0.01
ETAs = [1e-3,3.1e-3,1e-2,3.1e-2,1e-1,3.1e-1,1,3, 5.5, 10, 18, 30, 55]
T = 0.1

ff_list = np.zeros(M)
ff_list_e = np.zeros(M)

Es_th=4/3
Gamma_th = 6*SIZE/np.pi*np.sqrt(Es_th/2/np.pi/T)*np.exp(-Es_th/T)
def Langer(eta):
    o1 = np.sqrt(3)
    return (-eta/2+np.sqrt(eta**2/4+o1**2))/o1

#%%
""" Read data """

Times = []

path='Out/times/'
for i in range(M):
    eta = ETAs[i]
    if os.path.exists(path):
        files = os.listdir(path)
        timesList = []
        NumOfSim = 0
        for filename in files:
            if re.match('times_'+f'{DT:.6f}'+'_'+f'{T:.6f}'+'_'+f'{eta:.6f}',filename):
                filepath = os.path.join(path, filename)
                timesList += [np.loadtxt(filepath)]
        for j in range(0,len(timesList)): 
                NumOfSim += len(timesList[j])      
        times = np.zeros(NumOfSim)
        q = 0
        for l in range(0,len(timesList)): 
            for j in range(0,len(timesList[l])):
                times[q] = timesList[l][j]
                q += 1            
        Times += [times]

num_of_events = np.zeros(M)

for i in range(M):
    for j in range(len(Times[i])):
        if Times[i][j] == 0:
            Times[i][j] = TimeOfSim
    Times[i].sort()
    num_of_events[i] = (Times[i]<TimeOfSim-1).sum()

t_points = 10000
t_l = [None]*M
Psurv = np.zeros((M,t_points))
for j in range(M):        
    t_l[j] = np.linspace(0,TimeOfSim-0.1,t_points)
    for i in range(len(Times[j])):
        if Times[j][i] == 0:
            Times[j][i] = TimeOfSim
    Times[j].sort()
for j in range(M):
    for i in range(t_points):
        Psurv[j,i] = (Times[j]>t_l[j][i]).sum()/len(Times[j])

def Gamma1(i, tMin=500, log_P_max=2):
    log_Psurv = np.log(Psurv[i])   
    iMin = int(tMin*t_points/TimeOfSim)
    iMax = np.where((-log_Psurv<log_P_max))[0][-1]
    tMax = t_l[i][iMax]
    num = ((Times[i]>tMin)&(Times[i]<tMax)).sum()
    Gamma = - (np.log(Psurv[i,iMax])-np.log(Psurv[i,iMin]))/(tMax-tMin)
    Gamma_e = Gamma/np.sqrt(num)*np.sqrt((np.exp(Gamma*tMax)-1)/Gamma/tMax)
    return (Gamma, Gamma_e)
    
(Gamma, Gamma_e) = (np.zeros(M), np.zeros(M))
for i in range(M):
    (Gamma[i], Gamma_e[i]) = Gamma1(i)

ff_list[:] = Gamma[:]/Gamma_th
ff_list_e[:] = Gamma_e[:]/Gamma_th

#%%
""" Big plot """

t = np.linspace(1e-4, 2e2, 10000)
fig, ax = plt.subplots(figsize=(6,4))
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=16)
ax.axhline(y=1, linestyle='dashdot', linewidth=0.75, color='black')
ax.set_xlabel(r'$\hat{\eta}$',fontsize=18)
ax.set_ylabel(r'$A^{(sim)}/A_E$',fontsize=18)
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
plt.xscale('log')
plt.ylim((0, 1.05))
plt.xlim((7e-4,70))
ax.plot(t, Langer(t), 'red', linestyle='dashed', linewidth=1, dashes=[6,2])
ax.errorbar(eta_list, ff_list, ff_list_e, color='black', fmt='o', ms=5, markerfacecolor='black', elinewidth=1, capsize=2)
plt.savefig('Decay_rate_full.pdf', bbox_inches='tight')
plt.show()

#%%
""" Small plot """

def f(t,a,b):
    return a+t*b

ff_list_short = ff_list[7:] 
ff_list_e_short = ff_list_e[7:] 
eta_list_short = eta_list[7:]

params, cov = curve_fit(f,np.log10(eta_list_short[2:]),np.log10(ff_list_short[2:]),sigma=ff_list_e_short[2:]/ff_list_short[2:],absolute_sigma=True)
a = params[0]
b = params[1]

str1 = r'$y=$'+f'{b:.2f}'+r'$x+$'+f'{a:.2f}'
t = np.linspace(0.5*eta_list_short[0],2*eta_list_short[-1],1000)
ts = np.linspace(0.95*eta_list_short[2],2*eta_list_short[-1],1000)
fig, ax = plt.subplots(figsize=(6,4))
ax.set_xlabel(r'$\hat{\eta}$',fontsize=18)
ax.set_ylabel(r'$A^{(sim)}/A_E$',fontsize=18)
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=16)
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
plt.xscale('log',base=10)
plt.yscale('log',base=10)
ax.set_xticks([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80])
ax.set_xticklabels([None,None,None,None,5,None,None,None,None,10,20,None,None,50,None,None,None])
ax.errorbar(eta_list_short,ff_list_short,ff_list_e_short,color='black', fmt='o', ms=5, markerfacecolor='black', elinewidth=1, capsize=2)
ax.set_yticks([0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6])
ax.set_yticklabels([0.02,None,None,0.05,None,None,None,None,0.1,0.2,None,None,0.5,None])
ax.plot(t,Langer(t),color='red',linestyle='dashdot')
ax.plot(ts,ts**b*10**a,color='blue',linestyle='dashed',linewidth=1)
plt.xlim((0.8*eta_list_short[0],1.2*eta_list_short[-1]))
plt.ylim((1.5e-2,7e-1))
plt.text(2e1,2e-1,str1,fontsize=16,color='blue')
plt.savefig('Decay_rate_as.pdf', bbox_inches='tight')
plt.show()

#%%
