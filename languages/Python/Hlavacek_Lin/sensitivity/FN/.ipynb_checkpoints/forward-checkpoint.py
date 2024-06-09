import numpy as np
from scipy.integrate import solve_ivp
a = 0.7
b = 0.8
tau = 12.5
I = 1.3
tSpan = np.linspace(0, 200, 1000)
V = 0.0
W = 0.0
par = np.array([a,b,tau,I])
IC = state = np.array([V,W])
stateDim, parDim = len(IC), len(par)

def F(t,state,par):
    V,W = state
    a,b,tau,I = par
    return np.array([V - (V**3)/3.0 - W + I, (1/tau)*(V+a-b*W)])    
        
def J(t,state,par):
    V,W = state
    a,b,tau,I = par
    return np.array([[1-V**2,-1],[1/tau,-b/tau]])     
        
def dFdtheta_constant(t,state,par):
    V,W = state
    a,b,tau,I = par
    return np.array([[0,0,0,1],[1/tau,-W/tau,-(tau)**(-2)*(V+a-b*W),0]])
 
def jointF(t,jointState,par):

    x = jointState[:len(state)]
    s = jointState[len(state):].reshape((len(state),len(par)))
    dx = F(t,x,par)
    ds = (J(t,x,par).dot(s)+ dFdtheta_constant(t,x,par)).reshape((len(state)*len(par),))
    return np.hstack((dx,ds))
    
stateDim = len(state)
nn = len(par)
jointState = np.zeros(stateDim*(nn+1))
jointState[:stateDim] = IC[:]
sol_buffer = solve_ivp(fun=lambda t,z: jointF(t,z,par), t_span=(tSpan[0],tSpan[-1]),y0=jointState, t_eval=tSpan, method='LSODA', rtol=1e-6,atol=1e-6)
print(sol_buffer.y[2:6,-1])