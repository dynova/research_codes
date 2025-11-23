from numpy import *
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':12})
colors = cm.plasma(linspace(0,1,12))
parLog = load('/Users/amallela/Documents/bmab/inferences/Dallas4-n5-parLog.npy')
fig = plt.figure()
fig.set_size_inches((100,100))
ax = [plt.subplot(23,23,i+1) for i in range(529)]
onset = 220000
end = 340000
parLabels = ['${ts}_0$','${ts}_1$','${ts}_2$','${ts}_3$','${ts}_4$','${ts}_5$','$\\lambda_0$','$\\lambda_1$','$\\lambda_2$','$\\lambda_3$','$\\lambda_4$','$p_0$','$p_1$','$p_2$','$p_3$','$p_4$','$\\theta_0$','$\\theta_1$','$y_0$','$y_1$','$\\beta$','$f_D','$r$']
for i in range(23):
    for j in range(23):
        if i==j:
            yy,xx=histogram(log(parLog[onset:end,i]), bins=30)
            xc = 0.5*xx[:-1]+0.5*xx[1:]
            ax[i*23+j].bar(xc, yy/(xc[1]-xc[0])/sum(yy), width=xc[1]-xc[0], facecolor=colors[0],lw=0)
            ax[i*23+j].set_xlabel('$log('+parLabels[i][1:-1]+')$')
            ax[i*23+j].set_ylabel('Posterior')
            ax[i*23+j].ticklabel_format(axis='y', style='sci', scilimits=(-1,1),useMathText=True)
        else:
            data1 = log(parLog[onset:end,j])
            data2 = log(parLog[onset:end,i])
            ax[i*23+j].hist2d(data1,data2, cmap=cm.plasma,bins=30)
            ax[i*23+j].set_xlabel('$log('+parLabels[j][1:-1]+')$')
            ax[i*23+j].set_ylabel('$log('+parLabels[i][1:-1]+')$')
            ax[i*23+j].set_xticks([])
            ax[i*23+j].set_yticks([])
fig.tight_layout()
fig.savefig('Dallas-pairs-plot.pdf',dpi=1200)