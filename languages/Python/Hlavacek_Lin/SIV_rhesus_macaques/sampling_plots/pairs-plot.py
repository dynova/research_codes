from numpy import *
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':12})
colors = cm.viridis(linspace(0,1,20))
parLog = zeros((10000,20,18))

for ii in range(20):
    parLog[:,ii,:] = genfromtxt('/Users/amallela/Documents/perelson/sampling/output/Results/A_MCMC/Runs/params_'+str(ii)+'.txt')[33000:43000,:]    
parLog = parLog.reshape((200000,18))

fig = plt.figure()
fig.set_size_inches((100,100))
ax = [plt.subplot(18,18,i+1) for i in range(324)]
parLabels = ['$K_a A_0$','$K_x R_T$','$R_1$','$R_2$','$c$','$f$','$f_p$','$f_r$','$\\hat{N}_{c_{ref}}$','$\\hat{R}_{T_{ref}}$','$\\hat{a}_{ref}$','$\hat{\\mu}$','$\\hat{p}$','$k_E$','$k_f R_T$','${k}_{\\minus c}$','${k}_{\\plus c}$','$\\rho$']
for i in range(18):
    for j in range(18):
        if i==j:
            yy,xx=histogram(log10(parLog[:,i]), bins=50)
            xc = 0.5*xx[:-1]+0.5*xx[1:]
            ax[i*18+j].bar(xc, yy/(xc[1]-xc[0])/sum(yy), width=xc[1]-xc[0], facecolor=colors[0],lw=0)
            ax[i*18+j].set_xlabel(r'$\log_{10}('+parLabels[i][1:-1]+')$')
            ax[i*18+j].set_ylabel('Frequency')
            ax[i*18+j].ticklabel_format(axis='y', style='sci', scilimits=(-1,1),useMathText=True)
        else:
            data1 = log10(parLog[:,j])
            data2 = log10(parLog[:,i])
            ax[i*18+j].hist2d(data1,data2, cmap=cm.viridis,bins=50)
            ax[i*18+j].set_xlabel(r'$\log_{10}('+parLabels[j][1:-1]+')$')
            ax[i*18+j].set_ylabel(r'$\log_{10}('+parLabels[i][1:-1]+')$')
            ax[i*18+j].set_xticks([])
            ax[i*18+j].set_yticks([])
fig.tight_layout()
fig.savefig('pairs-plot.pdf',dpi=1200)