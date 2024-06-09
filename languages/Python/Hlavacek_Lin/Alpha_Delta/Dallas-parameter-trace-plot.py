from numpy import *
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':12})
colors = cm.plasma(linspace(0,1,12))
parLog = load('/Users/amallela/Documents/bmab/inferences/Dallas4-n5-parLog.npy')
fig = plt.figure()
fig.set_size_inches((100,100))
ax = [plt.subplot(6,4,i+1) for i in range(24)]
onset = 220000
end = 340000
parLabels = ['${ts}_0$','${ts}_1$','${ts}_2$','${ts}_3$','${ts}_4$','${ts}_5$','$\\lambda_0$','$\\lambda_1$','$\\lambda_2$','$\\lambda_3$','$\\lambda_4$','$p_0$','$p_1$','$p_2$','$p_3$','$p_4$','$\\theta_0$','$\\theta_1$','$y_0$','$y_1$','$\\beta$','$f_D','$r$']
for i in range(23):
	data = log(parLog[onset:end,i])
	ax[i].plot(data)
	ax[i].set_xlabel('Epoch')
	ax[i].set_ylabel('$log('+parLabels[i][1:-1]+')$')
fig.tight_layout()
fig.savefig('Dallas-parameter-trace-plot.pdf',dpi=1200)