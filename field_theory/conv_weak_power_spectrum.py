###########################
#  Weak convergence test  #
###########################
#
# The script reads field power spectra of stochastic solutions computed with different time steps.
# At each time step, it computes the average and the variance of the observable constructed from the long Fourier modes of the field.
# The linear fit is made in the log-log scale and the rate of convergence to the theory value is measured.

#%%
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
Ts = np.array([10.0, 9.8, 10.0, 10.0, 9.9, 10.0, 10.0, 9.8])

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
# read data

PS = []
PSE = []

GammaList = [2,3]

for w in range(2):
    
    GAMMA = GammaList[w]
    PS_av = np.zeros(N)
    PS_var = np.zeros(N)

    for i in range(N):
        
        ps_phi = []
        TIME_SPAN = Ts[i]
        N_min = int(N_SAMPLE*20*ETA/TIME_SPAN) # give it some time to thermalize        
        DT = DTs[i]
        str = 'ps_phi_'+f'{ETA:.6f}'+'_'+f'{TEMP:.6f}'+'_'+f'{DT:.6f}'+'_'+f'{BETA}'+'_'+f'{GAMMA}'+'_'+f'{SIGN}'

        for filename in files:
            if re.match(str,filename):
                filepath = os.path.join(path,filename)
                ps = np.loadtxt(filepath,usecols=range(k_max))

                for j in np.arange(N_min,N_SAMPLE,1):
                    ps_phi += [DX*np.sum(ps[j][:])/k_max]

        SimNo = len(ps_phi)
        ps_phi = np.array(ps_phi)
        PS_av[i] = ps_phi.sum()/SimNo
        for j in range(SimNo):
            PS_var[i] += (ps_phi[j]-PS_av[i])**2
        PS_var[i] = np.sqrt(PS_var[i])/(SimNo-1)

    PS += PS_av
    PSE += PS_var
  
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
ax.set_ylabel(r'$|\bar{{\cal S}}/{\cal S}_3 - 1|$',fontsize=18)
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

    params, cov = curve_fit(f, np.log10(DTs[-5:-1]), np.log10(np.abs(PS[w][-5:-1]/psTh[3]-1)))
    (a, b) = (params[0], params[1])
    print(a,b)

    if b>0:
        str1 = r'$y=$'+f'{a:.2f}'+r'$x+$'+f'{np.abs(b):.2f}'
    else:
        str1 = r'$y=$'+f'{a:.2f}'+r'$x-$'+f'{np.abs(b):.2f}'

    ax.plot(h, h**a*10**b,color=ColorList[w],linestyle=LineStyleList[w])
    ax.errorbar(DTs*(1+0.01*w), np.abs(PS[w]/psTh[3]-1), PSE[w]/psTh[3], color=ColorList1[w], fmt=markerList[w], ms=5, markerfacecolor=ColorList1[w], elinewidth=1, capsize=2)
    plt.text(xPos[w],yPos[w],str1,fontsize=16,color=ColorList[w])
plt.savefig('PS_conv.pdf', bbox_inches='tight')
#plt.show()

#%%
# accuracy plot

fig, ax = plt.subplots(figsize=(6,4))

ax.set_xlabel(r'$h$',fontsize=18)
ax.set_ylabel(r'$\bar{{\cal S}}/{\cal S}_0 - 1$',fontsize=18)
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

e1 = ax.errorbar(DTs, PS[0]/psTh[0]-1, PSE[0]/psTh[0], label=r'$(2,4)$', color='blue', fmt='o', ms=5, markerfacecolor='blue', elinewidth=1, capsize=2)
e2 = ax.errorbar(DTs*1.02, PS[1]/psTh[0]-1, PSE[1]/psTh[0], label=r'$(3,4)$', color='red', fmt='s', ms=5, markerfacecolor='red', elinewidth=1, capsize=2)
plt.text(0.53,-0.093,r'$\textit{1-loop}$',fontsize=14)
plt.text(0.53,-0.084,r'$\textit{2-loop}$',fontsize=14)
plt.text(0.53,-0.0868,r'$\textit{3-loop}$',fontsize=14)

ax.legend(handles=[e1, e2], loc='lower left',title=r'$(\gamma,s)$',alignment='left',fontsize=14,title_fontsize=14,markerfirst=False)
plt.savefig('PS_acc.pdf', bbox_inches='tight')
#plt.show()

# %%
