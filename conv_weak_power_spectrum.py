###########################
#  Weak convergence test  #
###########################
#
# The script reads field power spectra of stochastic solutions computed with different time steps.
# At each time step, it computes the average and the variance of the observable constructed from the long Fourier modes of the field.
# The linear fit is made in the log-log scale and the rate of convergence to the theory value is measured.

import re 
import os
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

PS = np.zeros(N)
PSE = np.zeros(N)

path = '/home/olborium/scratch/power_spectrum_phi/'
files = os.listdir(path)

for i in range(N):

    PSS = []

    DT = DTs[i]
    TIME_SPAN = DT*N_SAMPLE
    N_min = int(N_SAMPLE*10/ETA/TIME_SPAN) # give it some time to thermalize

    str = 'ps_'+f'{ETA:.6f}'+'_'+f'{TEMP:.6f}'+'_'+f'{DT:.6f}'+'_'+f'{BETA}'+'_'+f'{GAMMA}'+'_'+f'{SIGN}'
    for filename in files:
        if re.match(str,filename):
            filepath = os.path.join(path,filename)
            ps = np.loadtxt(filepath)

            for j in np.arange(N_min,N_SAMPLE,1):
                PSS += [DX*TEMP*np.sum(ps[j][0:k_max])/k_max]         
                                
    SimNo = len(PSS)
    PSS = np.array(PSS)
    PS[i] = PSS.sum()/SimNo

    for j in range(SimNo):
        PSE[i] += (PS[i] - PSS[j])**2
    PSE[i] = np.sqrt(PSE[i])/(SimNo-1)   

# make a fit:

def f(h,a,b):
    return a*h + b

params, cov = curve_fit(f, np.log10(DTs[-5:-2]), np.log10(PS[-5:-2]))
(a, b) = (params[0], params[1])
print(a,b)

# make a plot:

str =  '_'+f'{ETA:.3f}'+\
        '_'+f'{TEMP:.2f}'+\
        '_'+f'{BETA}'+\
        '_'+f'{GAMMA}'
if b>0:
    str1 = r'$y=$'+f'{a:.2f}'+r'$x+$'+f'{np.abs(b):.2f}'
else:
    str1 = r'$y=$'+f'{a:.2f}'+r'$x-$'+f'{np.abs(b):.2f}'
h = np.linspace(DTs[0],DTs[-1],100)
fig, ax = plt.subplots(figsize=(6,4))
ax.set_xlabel(r'$h$',fontsize=18)
ax.set_ylabel(r'$\Delta_k$',fontsize=18)
ax.set_xscale('log',base=10)
ax.set_yscale('log',base=10)
plt.ylim((1e-3,1e-1))
plt.xlim((6e-2,8e-1))
ax.set_xticks([0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
ax.set_xticklabels([None, None, None, 0.1, 0.2, None, 0.4, None, 0.6, None])
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=16)
ax.plot(h, h**a*10**b,'blue',linestyle='dashed')
ax.errorbar(DTs, PS, PSE, color='darkblue', fmt='o', ms=5, markerfacecolor='black', elinewidth=1, capsize=2)
plt.text(2e-1,3e-2,str1,fontsize=14,color='blue')
plt.savefig('PS'+str+'.pdf', bbox_inches='tight')
plt.clf()
