from numpy import *
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':12})
colors = cm.plasma(linspace(0,1,12))
parLog = load('/Users/amallela/Documents/bmab/inferences/Houston2-n4-parLog.npy')
fig = plt.figure()
fig.set_size_inches((100,100))
ax = [plt.subplot(5,4,i+1) for i in range(20)]
onset = 200000
end = 380000
parLabels = ['${ts}_0$','${ts}_1$','${ts}_2$','${ts}_3$','${ts}_4$','$\\lambda_0$','$\\lambda_1$','$\\lambda_2$','$\\lambda_3$','$p_0$','$p_1$','$p_2$','$p_3$','$\\theta_0$','$\\theta_1$','$y_0$','$y_1$','$\\beta$','$f_D','$r$']
for i in range(20):
	data = log(parLog[onset:end,i])
	ax[i].plot(data)
	ax[i].set_xlabel('Epoch')
	ax[i].set_ylabel('$log('+parLabels[i][1:-1]+')$')
fig.tight_layout()
fig.savefig('Houston-parameter-trace-plot.pdf',dpi=1200)