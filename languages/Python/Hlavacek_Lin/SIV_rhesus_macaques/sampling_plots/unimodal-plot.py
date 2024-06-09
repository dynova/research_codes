from numpy import *
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
import arviz as az
from matplotlib.ticker import FuncFormatter
plt.rcParams.update({'font.size':200})
colors = cm.viridis(linspace(0,1,20))
parLog = zeros((10000,20,18))

for ii in range(20):
    parLog[:,ii,:] = genfromtxt('/Users/amallela/Documents/perelson/sampling/output/Results/A_MCMC/Runs/params_'+str(ii)+'.txt')[33000:43000,:]    
parLog = parLog.reshape((200000,18))

arr = zeros((25,3))

arr[:18,0] = median(parLog,axis=0)
arr[18,0] = 1./arr[4,0]
arr[19,0] = arr[8,0]*1e12
arr[20,0] = arr[9,0]*1e12
arr[21,0] = arr[10,0]*1e6
arr[22,0] = 1./arr[11,0]*365.25/12.0
arr[23,0] = 1./arr[12,0]/86400.0
arr[24,0] = 1./arr[13,0]*3.0

arr[:18,1:] = az.hdi(parLog,hdi_prob=0.95)
for i in range(1,3):
    arr[18,i] = 1./arr[4,i]
    arr[19,i] = arr[8,i]*1e12    
    arr[20,i] = arr[9,i]*1e12    
    arr[21,i] = arr[10,i]*1e6
    arr[22,i] = 1./arr[11,i]*365.25/12.0
    arr[23,i] = 1./arr[12,i]/86400.0
    arr[24,i] = 1./arr[13,i]*3.0

fig = plt.figure()
fig.set_size_inches((100,100))
parLabels = ['$K_aA_0$ (composite body weight-independent parameter; dimensionless)','$K_xR_T$ (single-site equilibrium crosslinking constant for virion-receptor interaction multiplied by the total receptor number; dimensionless)','$R_1$ (dimensionless)','$R_2$ (dimensionless)', '$c$ (rate constant for clearance of free virus in the interstitial fluid of lymphoid tissue; d$^{-1}$)','$f$ (fraction determining the relative lifetimes of short- and long-lived productively infected cells; d$^{-1}$)','$f_p$ (fraction that determines the values of $\\alpha$ and $\\gamma$; dimensionless)', '$f_r$ (fraction that determines the value of $r$; dimensionless)','$\hat{N}_{c_{ref}}$ (population (in trillions) of CD4- mononuclear cells at reference body weight)','$\hat{R}_{T_{ref}}$ (number (in trillions) of receptors for reference body weight)','$\hat{a}_{ref}$ (total FDC surface area for reference body weight; m$^{2}$)','$\hat{\mu}$ (rate constant for death of uninfected target cells; month$^{-1}$)','$\hat{p}$ (rate constant characterizing generation of virus by productively infected cells; d$^{-1}$)','$k_E$ (rate constant for progression of exposed cells through stages of the eclipse phase of infection; d$^{-1}$)','$k_fR_T$ (composite body weight-independent parameter; d$^{-1}$)','$k_{\minus c}$ (rate constant for reversal of immobilization/trapping/isolation; d$^{-1}$)','$k_{\plus c}$ (rate constant for immobilization/trapping/isolation of virion-crosslinked receptor aggregates; d$^{-1}$)','$\\rho$ (ratio in steady-state of persistent viral infection; dimensionless)']
for i in range(18):
    if i == 12:   
        ax = [plt.subplot(1,1,1)]
        data = parLog[:,i]*86400
        _, bins = histogram(log10(data), bins=50)
        ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
        ax[0].set_xscale("log")    
        ax[0].set_xlabel(parLabels[i],fontsize=200,wrap=True)
        ax[0].set_ylabel('Frequency')
        ax[0].tick_params(length=60,width=10)
        ax[0].tick_params(which='minor', length=30,width=5)
        ax[0].yaxis.set_tick_params(which='minor', bottom=False)
        ax[0].axvline(x=arr[i,0]*86400, ls='-', color='r', lw=50)    
        ax[0].axvline(x=arr[i,1]*86400, ls=':', color='y', lw=50)    
        ax[0].axvline(x=arr[i,2]*86400, ls=':', color='y', lw=50)    
        for axis in [ax[0].xaxis, ax[0].yaxis]:
            formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
            axis.set_major_formatter(formatter)    
        fig.tight_layout()
        fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
        fig.clf()
    else:
        ax = [plt.subplot(1,1,1)]
        data = parLog[:,i]
        _, bins = histogram(log10(data), bins=50)
        ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
        ax[0].set_xscale("log")    
        ax[0].set_xlabel(parLabels[i],fontsize=200,wrap=True)
        ax[0].set_ylabel('Frequency')
        ax[0].tick_params(length=60,width=10)
        ax[0].tick_params(which='minor', length=30,width=5)
        ax[0].yaxis.set_tick_params(which='minor', bottom=False)
        ax[0].axvline(x=arr[i,0], ls='-', color='r', lw=50)    
        ax[0].axvline(x=arr[i,1], ls=':', color='y', lw=50)    
        ax[0].axvline(x=arr[i,2], ls=':', color='y', lw=50)    
        for axis in [ax[0].xaxis, ax[0].yaxis]:
            formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
            axis.set_major_formatter(formatter)    
        fig.tight_layout()
        fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
        fig.clf()            
    
parLabels = ['$1/c$ (lifetime for clearance of free virus in the interstitial fluid of lymphoid tissue; d)','$N_{c_{ref}}$ (population of CD4- mononuclear cells at reference body weight)','$R_{T_{ref}}$ (number of receptors for reference body weight)','$a_{ref}$ (total FDC surface area for reference body weight; mm$^{2}$)','$1/\mu$ (lifetime for death of uninfected target cells; d)','$1/p$ (lifetime for generation of virus by productively infected cells; d)','$\\tau_E$ (expected duration for progression of exposed cells through one stage of the eclipse phase of infection; d)']

i = 18
ax = [plt.subplot(1,1,1)]
data = 1./parLog[:,4]
_, bins = histogram(log10(data), bins=50)
ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
ax[0].set_xscale("log")    
ax[0].set_xlabel(parLabels[0],fontsize=200,wrap=True)
ax[0].set_ylabel('Frequency')
ax[0].tick_params(length=60,width=10)
ax[0].tick_params(which='minor', length=30,width=5)
ax[0].yaxis.set_tick_params(which='minor', bottom=False)
ax[0].axvline(x=arr[i,0], ls='-', color='r', lw=50)    
ax[0].axvline(x=arr[i,1], ls=':', color='y', lw=50)    
ax[0].axvline(x=arr[i,2], ls=':', color='y', lw=50)    
for axis in [ax[0].xaxis, ax[0].yaxis]:
    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
    axis.set_major_formatter(formatter)    
fig.tight_layout()
fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
fig.clf()    

i = 19
ax = [plt.subplot(1,1,1)]
data = parLog[:,8]*1e12
_, bins = histogram(log10(data), bins=50)
ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
ax[0].set_xscale("log")    
ax[0].set_xlabel(parLabels[1],fontsize=200,wrap=True)
ax[0].set_ylabel('Frequency')
ax[0].tick_params(length=60,width=10)
ax[0].tick_params(which='minor', length=30,width=5)
ax[0].yaxis.set_tick_params(which='minor', bottom=False)
ax[0].axvline(x=arr[i,0], ls='-', color='r', lw=50)    
ax[0].axvline(x=arr[i,1], ls=':', color='y', lw=50)       
ax[0].axvline(x=arr[i,2], ls=':', color='y', lw=50)    
fig.tight_layout()
fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
fig.clf() 

i = 20
ax = [plt.subplot(1,1,1)]
data = parLog[:,9]*1e12
_, bins = histogram(log10(data), bins=50)
ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
ax[0].set_xscale("log")    
ax[0].set_xlabel(parLabels[2],fontsize=200,wrap=True)
ax[0].set_ylabel('Frequency')
ax[0].tick_params(length=60,width=10)
ax[0].tick_params(which='minor', length=30,width=5)
ax[0].yaxis.set_tick_params(which='minor', bottom=False)
ax[0].axvline(x=arr[i,0], ls='-', color='r', lw=50)        
ax[0].axvline(x=arr[i,1], ls=':', color='y', lw=50)        
ax[0].axvline(x=arr[i,2], ls=':', color='y', lw=50)    
fig.tight_layout()
fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
fig.clf() 

i = 21
ax = [plt.subplot(1,1,1)]
data = parLog[:,10]*1e6
_, bins = histogram(log10(data), bins=50)
ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
ax[0].set_xscale("log")    
ax[0].set_xlabel(parLabels[3],fontsize=200,wrap=True)
ax[0].set_ylabel('Frequency')
ax[0].tick_params(length=60,width=10)
ax[0].tick_params(which='minor', length=30,width=5)
ax[0].yaxis.set_tick_params(which='minor', bottom=False)
ax[0].axvline(x=arr[i,0], ls='-', color='r', lw=50)    
ax[0].axvline(x=arr[i,1], ls=':', color='y', lw=50)        
ax[0].axvline(x=arr[i,2], ls=':', color='y', lw=50)       
fig.tight_layout()
fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
fig.clf() 

i = 22
ax = [plt.subplot(1,1,1)]
data = 1./parLog[:,11]*365.25/12.0
_, bins = histogram(log10(data), bins=50)
ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
ax[0].set_xscale("log")    
ax[0].set_xlabel(parLabels[4],fontsize=200,wrap=True)
ax[0].set_ylabel('Frequency')
ax[0].tick_params(length=60,width=10)
ax[0].tick_params(which='minor', length=30,width=5)
ax[0].yaxis.set_tick_params(which='minor', bottom=False)
ax[0].axvline(x=arr[i,0], ls='-', color='r', lw=50)    
ax[0].axvline(x=arr[i,1], ls=':', color='y', lw=50)    
ax[0].axvline(x=arr[i,2], ls=':', color='y', lw=50)  
for axis in [ax[0].xaxis, ax[0].yaxis]:
    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
    axis.set_major_formatter(formatter)    
fig.tight_layout()
fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
fig.clf() 

i = 23
ax = [plt.subplot(1,1,1)]
data = 1./parLog[:,12]/86400.0
_, bins = histogram(log10(data), bins=50)
ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
ax[0].set_xscale("log")    
ax[0].set_xlabel(parLabels[5],fontsize=200,wrap=True)
ax[0].set_ylabel('Frequency')
ax[0].tick_params(length=60,width=10)
ax[0].tick_params(which='minor', length=30,width=5)
ax[0].yaxis.set_tick_params(which='minor', bottom=False)
ax[0].axvline(x=arr[i,0], ls='-', color='r', lw=50)    
ax[0].axvline(x=arr[i,1], ls=':', color='y', lw=50)    
ax[0].axvline(x=arr[i,2], ls=':', color='y', lw=50) 
for axis in [ax[0].xaxis, ax[0].yaxis]:
    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
    axis.set_major_formatter(formatter)    
fig.tight_layout()
fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
fig.clf()

i = 24
ax = [plt.subplot(1,1,1)]
data = 1./parLog[:,13]*3.0
_, bins = histogram(log10(data), bins=50)
ax[0].hist(data, bins=10**bins, facecolor=colors[0],lw=0)
ax[0].set_xscale("log")    
ax[0].set_xlabel(parLabels[6],fontsize=200,wrap=True)
ax[0].set_ylabel('Frequency')
ax[0].tick_params(length=60,width=10)
ax[0].tick_params(which='minor', length=30,width=5)
ax[0].yaxis.set_tick_params(which='minor', bottom=False)
ax[0].axvline(x=arr[i,0], ls='-', color='r', lw=50) 
ax[0].axvline(x=arr[i,1], ls=':', color='y', lw=50)  
ax[0].axvline(x=arr[i,2], ls=':', color='y', lw=50)
for axis in [ax[0].xaxis, ax[0].yaxis]:
    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
    axis.set_major_formatter(formatter)    
fig.tight_layout()
fig.savefig('par'+str(i+1)+'.pdf',dpi=1200)
fig.clf()