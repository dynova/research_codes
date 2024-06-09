from numpy import *
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':200})
colors = cm.viridis(linspace(0,1,20))
fig = plt.figure()
fig.set_size_inches((100,100))
onset = 58000
end = 68000
parLabels = ['KaA0','KxRT','R1','R2','c','f','fp','f','hatNcref','hatRTref','hataref','hatmu','hatp','kE','kfRT','kmc','kpc','rho']
ax = [plt.subplot(1,1,1)]
for ii in range(20):
    logLLog = genfromtxt('/Users/amallela/Documents/perelson/sampling/output/Results/A_MCMC/Runs/scores_'+str(ii)+'.txt')    
    data = logLLog[onset:end]
    ax[0].plot(data,label='Chain '+str(ii+1),color=colors[ii])
    ax[0].set_xlabel('Epoch')
    ax[0].set_ylabel('Log-likelihood')
    ax[0].tick_params(length=60,width=10)
fig.tight_layout()
fig.savefig('markov-chain-trace-plot.pdf',dpi=1200)        