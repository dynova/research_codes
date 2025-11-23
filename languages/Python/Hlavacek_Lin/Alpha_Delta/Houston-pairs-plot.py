from numpy import *
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':12})
colors = cm.plasma(linspace(0,1,12))
parLog = load('/Users/amallela/Documents/bmab/inferences/Houston2-n4-parLog.npy')
fig = plt.figure()
fig.set_size_inches((100,100))
ax = [plt.subplot(20,20,i+1) for i in range(400)]
onset = 200000
end = 380000
parLabels = ['${ts}_0$','${ts}_1$','${ts}_2$','${ts}_3$','${ts}_4$','$\\lambda_0$','$\\lambda_1$','$\\lambda_2$','$\\lambda_3$','$p_0$','$p_1$','$p_2$','$p_3$','$\\theta_0$','$\\theta_1$','$y_0$','$y_1$','$\\beta$','$f_D','$r$']
for i in range(20):
    for j in range(20):
        if i==j:
            yy,xx=histogram(log(parLog[onset:end,i]), bins=30)
            xc = 0.5*xx[:-1]+0.5*xx[1:]
            ax[i*20+j].bar(xc, yy/(xc[1]-xc[0])/sum(yy), width=xc[1]-xc[0], facecolor=colors[0],lw=0)
            ax[i*20+j].set_xlabel('$log('+parLabels[i][1:-1]+')$')
            ax[i*20+j].set_ylabel('Posterior')
            ax[i*20+j].ticklabel_format(axis='y', style='sci', scilimits=(-1,1),useMathText=True)
        else:
            data1 = log(parLog[onset:end,j])
            data2 = log(parLog[onset:end,i])
            ax[i*20+j].hist2d(data1,data2, cmap=cm.plasma,bins=30)
            ax[i*20+j].set_xlabel('$log('+parLabels[j][1:-1]+')$')
            ax[i*20+j].set_ylabel('$log('+parLabels[i][1:-1]+')$')
            ax[i*20+j].set_xticks([])
            ax[i*20+j].set_yticks([])
fig.tight_layout()
fig.savefig('Houston-pairs-plot.pdf',dpi=1200)