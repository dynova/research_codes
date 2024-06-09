import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, ticker
n1, n2 = 100, 100
arr_eps = np.logspace(-4,0,n1)
arr_tol = np.logspace(-8,-4,n2)
eps, tol = np.meshgrid(arr_eps, arr_tol)
sens = np.load('run_contour.npy')
fig, ax = plt.subplots()
cs = ax.contourf(eps, tol, np.abs(sens), locator=ticker.LogLocator(), cmap=cm.plasma)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$\\epsilon$')
ax.set_ylabel('relative tolerance = absolute tolerance')
cbar = fig.colorbar(cs)
cbar.set_label('Sensitivity (magnitude)')
ax.set_title('Diffrax: Dopri8 with automatic differentiation')
plt.savefig("contour_plot.pdf")