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

SIZE = 100
N_x = 8192
N_SAMPLE = 52
BETA = 4

TEMP = 0.1
ETA = 1
SIGN = 1
N = 8
DTs = np.array([0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7])
Ts = np.array([0.5, 0.49, 0.5, 0.6, 0.6, 0.8, 1.0, 1.4])

DX = SIZE/N_x

k = np.linspace(0, np.pi/DX, N_x//2)

def f(x, a):
    return a
#%%

path = 'Out/power_spectrum_av/'
files = os.listdir(path)
j_trunc = 31

# theoretical power spectrum up to 3 loops:
def PSTh(k,n):
    return 1/(k**2+1 \
              + (n>0)*( 3*TEMP/2 \
              + (n>1)*(-9*TEMP**2/8-9*TEMP**2/2/(k**2+9) \
              + (n>2)*( 315*TEMP**3/128+27*(5*k**2+117)/16/(k**2+9)**2*TEMP**3 ) ) ) )

psTh = np.zeros(4)
for i in range(4):
    for j in range(j_trunc):
        kj = 2*np.pi*j/SIZE
        psTh[i] += PSTh(kj,i)   
    psTh[i] /= j_trunc

def f(h,a,b):
    return a*h + b

#%%
PS = []
PSE = []

GammaList = [2,3]

for w in range(2):
    
    ps_phi = []
    GAMMA = GammaList[w]

    for i in range(N):
        DT = DTs[i]
        str = 'ps_phi_'+f'{ETA:.6f}'+'_'+f'{TEMP:.6f}'+'_'+f'{DT:.6f}'+'_'+f'{BETA}'+'_'+f'{GAMMA}'+'_'+f'{SIGN}'+'_'
        pattern = re.compile(rf'^{str}\d+$')
        for filename in files:
            if pattern.search(filename):
                print('loading file ', filename)
                filepath = os.path.join(path,filename)
                ps_phi  += [np.loadtxt(filepath)]

    DeltaPS = np.zeros(N)
    SigmaDeltaPS = np.zeros(N)

    for i in range(N):

        TIME_SPAN = Ts[i]*N_SAMPLE
        DT = DTs[i]
        N_min = int(N_SAMPLE*20*ETA/TIME_SPAN) # initial thermalization time excluded

        """ Power spectrum """

        ps_sum = np.zeros((N_SAMPLE))
        for l in range(N_SAMPLE):
            for j in range(j_trunc):
                kj = 2*np.pi*j/SIZE
                ps_sum[l] += ps_phi[i][l][j]*DX/TEMP
            ps_sum[l] /= j_trunc

        ps_sum_av = 0
        for l in np.arange(N_min,N_SAMPLE,1):
            ps_sum_av += ps_sum[l]
        ps_sum_av /= (N_SAMPLE-N_min)

        ps_sum_var = 0
        for l in np.arange(N_min,N_SAMPLE,1):
            ps_sum_var += (ps_sum_av - ps_sum[l])**2
        ps_sum_var = np.sqrt(ps_sum_var)
        ps_sum_var /= (N_SAMPLE-N_min)
        
        DeltaPS[i] = ps_sum_av
        SigmaDeltaPS[i] = ps_sum_var

    PS += [DeltaPS]
    PSE += [SigmaDeltaPS]
    
#%%
# convergence plot

xPos = [0.2,0.35]
yPos = [1e-2,1.05e-3]
ColorList = ['blue', 'red']
ColorList1 = ['blue','red']
LineStyleList = ['dashed','dashdot']
markerList = ['o','s']

fig, ax = plt.subplots(figsize=(6,4))

ax.set_xlabel(r'$h$',fontsize=18)
ax.set_ylabel(r'$|\bar{S}/S_3 - 1|$',fontsize=18)
ax.set_xscale('log',base=10)
ax.set_yscale('log',base=10)
plt.ylim((2e-4,3e-2))
plt.xlim((6e-2,8e-1))
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=16)
ax.set_xticks([0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
ax.set_xticklabels([None, None, None, 0.1, 0.2, None, 0.4, None, 0.6, None])

h = np.linspace(DTs[0],DTs[-1],100)

for w in range(2):
    # (42)
    if w==0:
        params, cov = curve_fit(f, np.log10(DTs[-5:-1]), np.log10(np.abs(PS[w][-5:-1]/psTh[3]-1)))
        (a, b) = (params[0], params[1])
        print(a,b)

    # (43)
    if w==1:
        params, cov = curve_fit(f, np.log10(DTs[-4:-1]), np.log10(np.abs(PS[w][-4:-1]/psTh[3]-1)))
        (a, b) = (params[0], params[1])
        print(a,b)

    str =  '_'+f'{ETA:.3f}'+\
            '_'+f'{TEMP:.2f}'+\
            '_'+f'{BETA}'+\
            '_'+f'{GAMMA}'+\
            '_'+f'{SIGN}'
    if b>0:
        str1 = r'$y=$'+f'{a:.2f}'+r'$x+$'+f'{np.abs(b):.2f}'
    else:
        str1 = r'$y=$'+f'{a:.2f}'+r'$x-$'+f'{np.abs(b):.2f}'

    ax.plot(h, h**a*10**b,color=ColorList[w],linestyle=LineStyleList[w])
    ax.errorbar(DTs*(1+0.01*w), np.abs(PS[w]/psTh[3]-1), PSE[w]/psTh[3], color=ColorList1[w], fmt=markerList[w], ms=5, markerfacecolor=ColorList1[w], elinewidth=1, capsize=2)
    plt.text(xPos[w],yPos[w],str1,fontsize=16,color=ColorList[w])
plt.savefig('PS_conv.pdf', bbox_inches='tight')
plt.show()

#%%
# accuracy plot

fig, ax = plt.subplots(figsize=(6,4))

ax.set_xlabel(r'$h$',fontsize=18)
ax.set_ylabel(r'$\bar{S}/S_0 - 1$',fontsize=18)
ax.set_xscale('log',base=10)
plt.ylim((-0.108,-0.078))
plt.xlim((6e-2,8e-1))
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=16)

plt.hlines(psTh[1]/psTh[0]-1,6e-2,8e-1,'black',linestyle='dashed')
plt.hlines(psTh[2]/psTh[0]-1,6e-2,8e-1,'black',linestyle='dashdot')
plt.hlines(psTh[3]/psTh[0]-1,6e-2,8e-1,'black',linestyle='dotted')

ax.set_xticks([0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
ax.set_xticklabels([None, None, None, 0.1, 0.2, None, 0.4, None, 0.6, None])

e1 = ax.errorbar(DTs, PS[0]/psTh[0]-1, PSE[0]/psTh[0], label=r'$(2,4)$', color='blue', fmt='s', ms=5, markerfacecolor='blue', elinewidth=1, capsize=2)
e2 = ax.errorbar(DTs*1.02, PS[1]/psTh[0]-1, PSE[1]/psTh[0], label=r'$(3,4)$', color='red', fmt='o', ms=5, markerfacecolor='red', elinewidth=1, capsize=2)

ax.legend(handles=[e1, e2], loc='lower left',title=r'$(\gamma,s)$',alignment='left',fontsize=14,title_fontsize=14,markerfirst=False)
plt.savefig('PS_acc.pdf', bbox_inches='tight')
plt.show()

# %%
