from numpy import *
from scipy.stats import loguniform
from scipy.integrate import solve_ivp
from scipy.special import loggamma
from scipy.optimize import minimize
from sys import argv
from numba import njit
from matplotlib import cm
import random
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':24})
colors = cm.plasma(linspace(0,1,12))
random.seed(10)
temp = load('/Users/amallela/Documents/bmab/inferences/Dallas4-n5-logLLog.npy')
fig = plt.figure(figsize=(24,24))
ax = plt.subplot(1,1,1)
ax.plot(temp[220000:340000])
ax.set_title('Dallas', fontsize = 30)
ax.set_xlabel('Epoch', fontsize = 24)
ax.set_ylabel('Log-likelihood')
fig.tight_layout()
plt.savefig('Dallas-markov-chain-trace-plot.pdf',dpi=1200)
plt.close()