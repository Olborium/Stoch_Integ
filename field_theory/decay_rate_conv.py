# Read the decay times, calculate the decay rate from the slope of the survival probability.
# Make figure 6 from the paper
#%%
import os
import re
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})
#%%
""" Parameters of the simulation """

SIZE = 100
N_x = 8192
DX = SIZE/N_x
TimeOfSim = 4000
M = 9
DTs = [0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.10]
eta = 10
T = 0.1

Es_th=4/3
Gamma_th = 6*SIZE/np.pi*np.sqrt(Es_th/2/np.pi/T)*np.exp(-Es_th/T)

o1 = np.sqrt(3)
Langer = (-eta/2+np.sqrt(eta**2/4+o1**2))/o1

#%%
""" read data """

pathList = ['Out/Gamma42/','Out/Gamma43/']
t_l_tot = []
Psurv_tot = []
Pref = np.zeros((2,M))
Pref_e = np.zeros((2,M))

for w in range(2):

  Times = []

  path=pathList[w]
  for i in range(M):
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

  def Gamma1(i, tMin=250, log_P_max=2):
      log_Psurv = np.log(Psurv[i])   
      iMin = int(tMin*t_points/TimeOfSim)
      iMax = np.where((-log_Psurv<log_P_max))[0][-1]
      tMax = t_l[i][iMax]
      num = ((Times[i]>tMin)&(Times[i]<tMax)).sum()
      Gamma = - (np.log(Psurv[i,iMax])-np.log(Psurv[i,iMin]))/(tMax-tMin)
      Gamma_e = Gamma/np.sqrt(num) # Simple error estimation
  #   Gamma_e = Gamma/np.sqrt(num)*np.sqrt((np.exp(Gamma*tMax)-1)/Gamma/tMax) # More precise error estimation, see Appendix D of 2408.06411
      return (Gamma, Gamma_e)
    
  (Gamma, Gamma_e) = (np.zeros(M), np.zeros(M))
  for i in range(M):
      (Gamma[i], Gamma_e[i]) = Gamma1(i)
  
  t_l_tot += [np.copy(t_l[2])]
  Psurv_tot += [np.copy(Psurv[2])]

  Pref[w,:] = Gamma[:]/Gamma_th
  Pref_e[w,:] = Gamma_e[:]/Gamma_th

#%%
# Decay rate as a function of h

t=np.linspace(0.015,0.12,1000)
fig, ax = plt.subplots(figsize=(6,4))
ax.set_xlabel(r'$h$',fontsize=18)
ax.set_ylabel(r'$\Gamma^{\rm{(sim)}}/\Gamma_{\rm{stat}}$',fontsize=18)
plt.xscale('log',base=10)
ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
ax.set_xticks([0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1])
ax.set_xticklabels([0.005,None,None,None,None,0.01,0.02,None,None,0.05,None,None,None,None,0.1])
e1 = ax.errorbar(DTs,Pref[0]/Langer,Pref_e[0]/Langer, label=r'$(2,4)$', color='blue', fmt='o', ms=5, markerfacecolor='blue', elinewidth=1, capsize=2)
e2 = ax.errorbar(DTs*1.02,Pref[1]/Langer,Pref_e[1]/Langer, label=r'$(3,4)$', color='red', fmt='s', ms=5, markerfacecolor='red', elinewidth=1, capsize=2)
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=14)
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.legend(handles=[e1, e2], loc='lower left',title=r'$(\gamma,s)$',alignment='left',fontsize=14,title_fontsize=14,markerfirst=False)
plt.ylim((0.3,1.0))
plt.savefig('Pref_conv.pdf', bbox_inches='tight')
plt.show()

#%%
# Probability plots

t = np.linspace(0, 4000, 1000)
fig, ax = plt.subplots(figsize=(6,4))
plt.ylim((-0.25,0))
plt.xlim((0,4000))  
t1 = ax.plot(t_l_tot[0], np.log(Psurv_tot[0]), 'blue',label=r'$(2,4)$')
t2 = ax.plot(t_l_tot[1], np.log(Psurv_tot[1]), 'red',linestyle='dashed',label=r'$(3,4)$')
ax.plot(t, -Langer*Gamma_th*t,'black',linestyle='dashdot' )
ax.tick_params(axis='both',which='both',right=True,top=True,direction='in',labelsize=14)
ax.set_xlabel(r'$t$',fontsize=18)
ax.set_ylabel(r'$\ln P_{\rm{surv}}$',fontsize=18)
ax.xaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.yaxis.grid(True,linestyle='dashed',linewidth=0.5,color='grey')
ax.legend(loc='lower left',title=r'$(\gamma,s)$',alignment='left',fontsize=14,title_fontsize=14,markerfirst=False)
plt.text(3100,-0.04,r'$h=0.01$',fontsize=16)
plt.text(2700,-0.227,r'$\textit{theory}$',fontsize=16,rotation=-39)
plt.savefig('Psurv.pdf', bbox_inches='tight')
plt.show()
#%%
