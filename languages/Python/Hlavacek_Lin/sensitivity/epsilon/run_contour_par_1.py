import numpy as np
import diffrax
import jax
import jax.numpy as jnp
import jax.scipy.stats as stats
from jax.scipy.special import gammaln
from time import time
from jax import grad
import matplotlib.pyplot as plt
jax.config.update("jax_enable_x64", True)
x = 1.
y = 1.
IC = state = jnp.array([x,y])
n1, n2 = 10, 10
temp = np.zeros((n1,n2))
arr_par1 = jnp.logspace(-4,0,n1)
eps = 1e-4
arr_tol = jnp.logspace(-8,-4,n2)
tf = 1.
tSpan = jnp.linspace(0, tf, 100)
cc = -1
@jax.jit
def RHS(t,state,par):
    x,y = state
    par1 = par
    return jnp.array([-(2+eps**(-1))*x+eps**(-1)*y**2,x-y*(y+par1)])
@jax.jit
def logL(par):
    terms = diffrax.ODETerm(RHS)
    solver = diffrax.Dopri8()
    t0 = tSpan[0]
    t1 = tSpan[-1]
    dt0 = None
    y0 = state
    saveat = diffrax.SaveAt(ts=tSpan)
    stepsize_controller = diffrax.PIDController(rtol=tol, atol=tol)       
    sol = diffrax.diffeqsolve(
    terms,
    solver,
    t0,
    t1,
    dt0,
    y0,
    args = par,
    saveat=saveat,
    stepsize_controller=stepsize_controller,
    max_steps = int(1e12),    
    )
    logLval = sol.ys[-1,-1]
    return logLval    
for ii in range(n1):
    par = arr_par1[ii]
    for jj in range(n2):
        tol = arr_tol[jj]
        cc += 1
        print(cc)
        temp[ii,jj] = grad(logL)(par)
np.save('run_contour_par_1.npy',temp)