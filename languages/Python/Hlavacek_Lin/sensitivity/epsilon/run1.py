### EULER ###
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
temp1 = np.zeros(9)
temp2 = np.zeros(9)
count = -1
for eps in jnp.logspace(-8,0,9):
    par = eps
    tf = 1.
    tSpan = jnp.linspace(0, tf, 100)
    def RHS(t,state,par):
        x,y = state
        eps = par
        return jnp.array([-(2+eps**(-1))*x+eps**(-1)*y**2,x-y*(y+1)])
    def logL(par):
        terms = diffrax.ODETerm(RHS)
        solver = diffrax.Euler()
        t0 = tSpan[0]
        t1 = tSpan[-1]
        dt0 = 1e-5
        y0 = state
        saveat = diffrax.SaveAt(ts=tSpan)
        sol = diffrax.diffeqsolve(
        terms,
        solver,
        t0,
        t1,
        dt0,
        y0,
        args = par,
        saveat=saveat,
        max_steps = int(1e12),    
        )
        logLval = sol.ys[-1,-1]
        return logLval    
    tolerance = 1e-6
    count += 1
    tic = time()
    temp1[count] = grad(logL)(par)
    temp2[count] = time()-tic
np.save('ad_euler_sens.npy',temp1)
np.save('ad_euler_time.npy',temp2)

### DOPRI8 ###
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
temp1 = np.zeros(9)
temp2 = np.zeros(9)
count = -1
for eps in jnp.logspace(-8,0,9):
    par = eps
    tf = 1.
    tSpan = jnp.linspace(0, tf, 100)
    def RHS(t,state,par):
        x,y = state
        eps = par
        return jnp.array([-(2+eps**(-1))*x+eps**(-1)*y**2,x-y*(y+1)])
    def logL(par):
        terms = diffrax.ODETerm(RHS)
        solver = diffrax.Dopri8()
        t0 = tSpan[0]
        t1 = tSpan[-1]
        dt0 = 1e-5
        y0 = state
        saveat = diffrax.SaveAt(ts=tSpan)
        stepsize_controller = diffrax.PIDController(rtol=tolerance, atol=tolerance)       
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
    tolerance = 1e-6
    count += 1
    tic = time()
    temp1[count] = grad(logL)(par)
    temp2[count] = time()-tic
np.save('ad_dopri8_sens.npy',temp1)
np.save('ad_dopri8_time.npy',temp2)
   
### FORWARD LSODA ###    
import numpy as np
from scipy.integrate import solve_ivp
temp1 = np.zeros(9)
temp2 = np.zeros(9)
count = -1
for eps in np.logspace(-8,0,9):
    tSpan = np.linspace(0, 1, 100)
    x = 1.
    y = 1.
    par = eps
    IC = state = np.array([x,y])
    stateDim, parDim = len(IC), 1
    def F(t,state,par):
        x,y = state
        eps = par
        return np.array([-(2+eps**(-1))*x+eps**(-1)*y**2,x-y*(y+1)])      
    def J(t,state,par):
        x,y = state
        eps = par
        return np.array([[-(2+eps**(-1)),2*eps**(-1)*y],[1,-2*y-1]])      
    def dFdtheta_constant(t,state,par):
        x,y = state
        eps = par
        return np.array([[eps**(-2)*x-eps**(-2)*y**2],[0]]) 
    def jointF(t,jointState,par):
        x = jointState[:len(state)]
        s = jointState[len(state):].reshape((len(state),1))
        dx = F(t,x,par)
        ds = (J(t,x,par).dot(s)+ dFdtheta_constant(t,x,par)).reshape((len(state),))
        return np.hstack((dx,ds))    
    stateDim = len(state)
    nn = 1
    jointState = np.zeros(stateDim*(nn+1))
    jointState[:stateDim] = IC[:]
    count += 1
    tic = time()
    sol_buffer = solve_ivp(fun=lambda t,z: jointF(t,z,par), t_span=(tSpan[0],tSpan[-1]),y0=jointState, t_eval=tSpan, method='LSODA', rtol=1e-6,atol=1e-6)
    temp1[count] = sol_buffer.y[-1,-1]
    temp2[count] = time()-tic
np.save('forward_lsoda_sens.npy',temp1)
np.save('forward_lsoda_time.npy',temp2)

### FORWARD Radau ###    
import numpy as np
from scipy.integrate import solve_ivp
temp1 = np.zeros(9)
temp2 = np.zeros(9)
count = -1
for eps in np.logspace(-8,0,9):
    tSpan = np.linspace(0, 1, 100)
    x = 1.
    y = 1.
    par = eps
    IC = state = np.array([x,y])
    stateDim, parDim = len(IC), 1
    def F(t,state,par):
        x,y = state
        eps = par
        return np.array([-(2+eps**(-1))*x+eps**(-1)*y**2,x-y*(y+1)])      
    def J(t,state,par):
        x,y = state
        eps = par
        return np.array([[-(2+eps**(-1)),2*eps**(-1)*y],[1,-2*y-1]])      
    def dFdtheta_constant(t,state,par):
        x,y = state
        eps = par
        return np.array([[eps**(-2)*x-eps**(-2)*y**2],[0]]) 
    def jointF(t,jointState,par):
        x = jointState[:len(state)]
        s = jointState[len(state):].reshape((len(state),1))
        dx = F(t,x,par)
        ds = (J(t,x,par).dot(s)+ dFdtheta_constant(t,x,par)).reshape((len(state),))
        return np.hstack((dx,ds))    
    stateDim = len(state)
    nn = 1
    jointState = np.zeros(stateDim*(nn+1))
    jointState[:stateDim] = IC[:]
    count += 1
    tic = time()
    sol_buffer = solve_ivp(fun=lambda t,z: jointF(t,z,par), t_span=(tSpan[0],tSpan[-1]),y0=jointState, t_eval=tSpan, method='Radau', rtol=1e-6,atol=1e-6)
    temp1[count] = sol_buffer.y[-1,-1]
    temp2[count] = time()-tic
np.save('forward_radau_sens.npy',temp1)
np.save('forward_radau_time.npy',temp2)