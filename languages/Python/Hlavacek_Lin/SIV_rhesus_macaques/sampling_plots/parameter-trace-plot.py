from numpy import *
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':48})
colors = cm.viridis(linspace(0,1,20))
fig = plt.figure()
fig.set_size_inches((100,100))
ax = [plt.subplot(6,3,jj+1) for jj in range(18)]
onset = 33000
end = 43000
parLabels = ['$K_a A_0$','$K_x R_T$','$R_1$','$R_2$','$c$','$f$','$f_p$','$f_r$','$\\hat{N}_{c_{ref}}$','$\\hat{R}_{T_{ref}}$','$\\hat{a}_{ref}$','$\hat{\\mu}$','$\\hat{p}$','$k_E$','$k_f R_T$','${k}_{\\minus c}$','${k}_{\\plus c}$','$\\rho$']
for ii in range(20):
    parLog = genfromtxt('/Users/amallela/Documents/perelson/sampling/output/Results/A_MCMC/Runs/params_'+str(ii)+'.txt')    
    for jj in range(18):
        data = log10(parLog[onset:end,jj])
        ax[jj].plot(data,color=colors[ii])
        ax[jj].set_xlabel('Epoch')
        ax[jj].set_ylabel(r'$\log_{10}('+parLabels[jj][1:-1]+')$')
        ax[jj].tick_params(length=30,width=5)
fig.tight_layout()
fig.savefig('parameter-trace-plot.pdf',dpi=1200)