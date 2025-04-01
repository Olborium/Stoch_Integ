###########################
#  Weak convergence test  #
###########################
#
# The script reads momenta power spectra of stochastic solutions computed with different time steps.
# At each time step, it computes the average and the variance of the effective temperature of long and short Fourier modes.
# The linear fit is made in the log-log scale and the rate of convergence to the theory value TEMP is measured.

import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})

SIZE = float(sys.argv[1])
N_x = int(sys.argv[2])
N_SAMPLE = int(sys.argv[3])
BETA = int(sys.argv[4])
GAMMA = int(sys.argv[5])
TEMP = float(sys.argv[6])
ETA = float(sys.argv[7])
N = int(sys.argv[8])
DTs = np.zeros(N)
Ts = np.zeros(N)
for i in range(N):
    DTs[i] = float(sys.argv[i+9])
    Ts[i] = float(sys.argv[i+9+N])
DX = SIZE/N_x

k = np.linspace(0, np.pi/DX, N_x//2)
k_max = int(SIZE/np.pi)

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
    N_min = int(N_SAMPLE*10/ETA/TIME_SPAN)  # give it some time to thermalize

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

# make a fit:

def f(h,a,b):
    return a*h + b

params, cov = curve_fit(f, np.log10(DTs[-4:]), np.log10(np.abs(TeffLon/TEMP-1)[-4:]))
(aL, bL) = (params[0], params[1])
print(aL,bL)

params, cov = curve_fit(f, np.log10(DTs[0:4]), np.log10(np.abs(TeffSho/TEMP-1)[0:4]))
(aS, bS) = (params[0], params[1])
print(aS,bS)

# make a plot:

str =  '_'+f'{ETA:.3f}'+\
        '_'+f'{TEMP:.2f}'+\
        '_'+f'{BETA}'+\
        '_'+f'{GAMMA}'

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
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=16)
ax.set_xlabel(r'$h$',fontsize=18)
ax.set_ylabel(r'$\delta_{\mathcal{T}} $',fontsize=18)
ax.set_xscale('log',base=10)
ax.set_yscale('log',base=10)
plt.ylim((1e-5,3e-0))
ax.plot(h, h**aS*10**bS,'red',linestyle='dashed')
ax.plot(h, h**aL*10**bL,'blue',linestyle='dashdot')
ax.errorbar(DTs, np.abs(TeffSho/TEMP-1), (TeffShoE/TEMP), color='darkred', fmt='o', ms=5, markerfacecolor='black', elinewidth=1, capsize=2)
ax.errorbar(DTs, np.abs(TeffLon/TEMP-1), (TeffLonE/TEMP), color='darkblue', fmt='D', ms=5, markerfacecolor='black', elinewidth=1, capsize=2)
plt.text(2.5e-2,2e-1,str2,fontsize=14,color='red')
plt.text(1e-1,2e-3,str1,fontsize=14,color='blue')
plt.savefig('Conv_weak_T'+str+'.pdf', bbox_inches='tight')
