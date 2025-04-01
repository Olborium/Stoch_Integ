###########################
#  Weak convergence test  #
###########################
#
# The script reads momenta power spectra of stochastic solutions computed with different time steps.
# At each time step, it computes the average and variance of the effective temperature of long and short Fourier modes.
# The linear fit is made in the log-log scale and the rate of convergence to the theory value TEMP is measured.

#%%
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

#%%
SIZE = 100
N_x = 8192
k_max = 31
N_SAMPLE = 50
BETA = 4
GAMMA = 3
TEMP = 0.1
ETA = 1
SIGN = 1

N = 14
DTs = np.array([0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7])
Ts = np.array([10.0, 9.8, 10.0, 10.0, 9.9, 10.0, 10.0, 9.8, 10.0, 10.09, 9.9, 10.0, 10.0, 9.8])

DX = SIZE/N_x

k = np.linspace(0, np.pi/DX, N_x//2)
#%%
TeffLon = np.zeros(N)
TeffLonE = np.zeros(N)
TeffSho = np.zeros(N)
TeffShoE = np.zeros(N)

path = '/home/olborium/scratch/power_spectrum_chi/'
files = os.listdir(path)

for i in range(N):

    TeffLonS = []
    TeffShoS = []

    DT = DTs[i]
    TIME_SPAN = DT*N_SAMPLE
    N_min = int(N_SAMPLE*20/ETA/TIME_SPAN)  # give it some time to thermalize

    str = 'ps_'+f'{ETA:.6f}'+'_'+f'{TEMP:.6f}'+'_'+f'{DT:.6f}'+'_'+f'{BETA}'+'_'+f'{GAMMA}'+'_'+f'{SIGN}'
    for filename in files:
        if re.match(str,filename):
            filepath = os.path.join(path,filename)
            ps = np.loadtxt(filepath)

            for j in np.arange(N_min,N_SAMPLE,1):
                TeffLonS += [DX*np.sum(ps[j][0:k_max])/k_max]        
                TeffShoS += [DX*np.sum(ps[j][N_x//2-k_max:N_x//2])/k_max]  
                                
    SimNo = len(TeffLonS)
    TeffLonS = np.array(TeffLonS)
    TeffShoS = np.array(TeffShoS)
    TeffLon[i] = TeffLonS.sum()/SimNo
    TeffSho[i] = TeffShoS.sum()/SimNo

    for j in range(SimNo):
        TeffLonE[i] += (TeffLon[i] - TeffLonS[j])**2
        TeffShoE[i] += (TeffSho[i] - TeffShoS[j])**2
    TeffLonE[i] = np.sqrt(TeffLonE[i])/(SimNo-1)
    TeffShoE[i] = np.sqrt(TeffShoE[i])/(SimNo-1)    

#%% 
# make a fit:

def f(h,a,b):
    return a*h + b

params, cov = curve_fit(f, np.log10(DTs[-4:]), np.log10(np.abs(TeffLon/TEMP-1)[-4:]))
(aL, bL) = (params[0], params[1])
print(aL,bL)

params, cov = curve_fit(f, np.log10(DTs[0:4]), np.log10(np.abs(TeffSho/TEMP-1)[0:4]))
(aS, bS) = (params[0], params[1])
print(aS,bS)

#%%
# make a plot:

if bL>0:
    str1 = r'$y=$'+f'{aL:.2f}'+r'$x+$'+f'{np.abs(bL):.2f}'
else:
    str1 = r'$y=$'+f'{aL:.2f}'+r'$x-$'+f'{np.abs(bL):.2f}'
if bS>0:
    str2 = r'$y=$'+f'{aS:.2f}'+r'$x+$'+f'{np.abs(bS):.2f}'
else:
    str2 = r'$y=$'+f'{aS:.2f}'+r'$x-$'+f'{np.abs(bS):.2f}'
h = np.linspace(DTs[0],DTs[-1],100)
fig, ax = plt.subplots(figsize=(6,4))

ax.set_xlabel(r'$h$',fontsize=18)
ax.set_ylabel(r'$\delta_{\mathcal{T}} $',fontsize=18)
ax.set_xscale('log',base=10)
ax.set_yscale('log',base=10)
plt.ylim((1e-5,3e-0))
ax.plot(h, h**aS*10**bS,'tomato',linestyle='dashed')
ax.plot(h, h**aL*10**bL,'blue',linestyle='dashdot')
ax.errorbar(DTs, np.abs(TeffSho/TEMP-1), TeffShoE/TEMP, color='red', fmt='o', ms=5, markerfacecolor='red', elinewidth=1, capsize=2)
ax.errorbar(DTs, np.abs(TeffLon/TEMP-1), TeffLonE/TEMP, color='blue', fmt='s', ms=5, markerfacecolor='blue', elinewidth=1, capsize=2)
plt.text(2.5e-2,2e-1,str2,fontsize=16,color='red')
plt.text(0.09,2e-3,str1,fontsize=16,color='blue')
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=16)
ax.set_yticks([1e-5,1e-4,1e-3,1e-2,1e-1,1e0])
ax.set_xticks([0.005,0.01,0.02,0.05,0.1,0.2,0.5])
ax.set_xticklabels([0.005,0.01,0.02,0.05,0.1,0.2,0.5])
plt.savefig('Temp43.pdf', bbox_inches='tight')
