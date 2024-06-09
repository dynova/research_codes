import numpy as np
from scipy.integrate import solve_ivp
np.random.seed(123)

nt = 101
tSpan = np.linspace(0, nt-1, nt)
tSpanSimulation = np.linspace(0, nt-1, nt)
ns = 50
ne = 70
N = 19216182.
t0 = 10.
beta1 = 1.
gamma = 0.2
f = 0.9
sigma = 0.1
par = np.array([t0, beta1, gamma, f, sigma])
s0 = 1 - 1/N
i0 = 1/N
r0 = 0
s = s0
i = i0
r = r0
IC = state = np.array([s,i,r])
stateDim, parDim = len(IC), len(par)

def F(t,state,par,ifOn=True):

    if ifOn:
        s, i, r = state
        t0, beta1, gamma, f, sigma = par    
        beta = beta1
        return np.array([-beta*s*i,beta*s*i-gamma*i,gamma*i])
    else:
        return np.zeros(np.shape(state))
        
def J(t,state,par,ifOn=True):

    if ifOn:
        s, i, r = state
        t0, beta1, gamma, f, sigma = par    
        beta = beta1
        return np.array([[-beta1*i,-beta1*s,0],[beta1*i,beta1*s-gamma,0],[0,gamma,0]])
    else:
        return np.zeros((len(state),len(state)))        
        
def dFdtheta_constant(t,state,par,ifOn=True):

    if ifOn:
        s, i, r = state
        t0, beta1, gamma, f, sigma = par    
        beta = beta1
        return np.array([[0,-s*i,0,0,0],[0,s*i,-i,0,0],[0,0,i,0,0]])
    else:     
        return np.zeros((len(state),len(par)))   
        
def d2Fdthetadx_constant(t,state,par,ifOn=True):

    if ifOn:
        s, i, r = state
        t0, beta1, gamma, f, sigma = par    
        beta = beta1
        return np.array([[0, 0, 0], [-i, -s, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [i, s, 0], [0, -1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 1, 0], [0, 0, 0], [0, 0, 0]])
    else:     
        return np.zeros((len(state)*len(par)),len(state))           
        
def jointF(t,jointState,par,ifOn=True):

    if ifOn:

        x = jointState[:len(state)]
        s = jointState[len(state):].reshape((len(state),len(par)))

        dx = F(t,x,par,ifOn)
        ds = (J(t,x,par,ifOn).dot(s)+ dFdtheta_constant(t,x,par,ifOn)).reshape((len(state)*len(par),))

        return np.hstack((dx,ds))

    else:

        return np.zeros(np.shape(jointState)) 
        
def jointJacobian(t,jointState,par,ifOn=True): 
    
    jointJ = np.zeros((stateDim+stateDim*nn,stateDim+stateDim*nn))
    
    if ifOn:
    
        x = jointState[:stateDim]
        s = jointState[stateDim:].reshape((stateDim,nn))           
    
        jointJ[:stateDim, :stateDim] = J(t,x,par,ifOn)
        jointJ[stateDim:, :stateDim] = d2Fdthetadx_constant(t,x,par,ifOn).reshape((stateDim*nn,stateDim))
        jointJ[stateDim:,stateDim:] = np.tensordot(J(t,x,par,ifOn), np.eye(nn), axes=0).swapaxes(1, 2).reshape((stateDim*nn,stateDim*nn))
    
        return jointJ
    
    else:
    
        return np.zeros((np.shape(jointState), np.shape(jointState)))                             

def process(par):
    sol = solve_ivp(fun=lambda t,z: F(t,z,par),jac=lambda t,z: J(t,z,par), t_span=(tSpan[0],tSpan[-1]),y0=state, t_eval=tSpan, method='LSODA',rtol=1e-6,atol=1e-6)
    newCases = par[-2]*sol.y[1,:]
    return newCases

def logL(par):
    output = process(par)
    sigma = par[-1]
    newCases = np.zeros(ne-ns)
    logLval = 0
    for i in range(ns,ne):
        newCases[i-ns] = output[i] + np.random.normal(output[i],sigma,1)
        logLval += -0.5*np.log(2*np.pi)-np.log(sigma)-0.5*(sigma**-2)*(output[i]-newCases[i-ns])**2
    return logLval, newCases

def deltaType3(state,par,ifOn=True):            
    s, i, r = state
    t0, beta1, gamma, f, sigma = par    

    if ifOn:
        beta = beta1
        return np.array([beta*s*i,-beta*s*i,0])
    else:
        return -F(t0,state,par,ifOn=True)
        
def forwardSens(par, parIndex):
    stateDim = len(state)
    nn = len(parIndex)
    jointState = np.zeros(stateDim*(nn+1))
    jointState[:stateDim] = IC[:]

    durations = np.array(par[:1])
    change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
    times = []

    for i in range(len(durations)+1):

        lower = change_times[i]
        upper = change_times[i+1]

        times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper)))


    for i in range(len(times)):
        if i==1:
            x,s = jointState[:stateDim],jointState[stateDim:].reshape((stateDim, nn))
            buffer = deltaType3(x, par, ifOn=True)
            s[:,0] += buffer    
            jointState = np.hstack((x, s.reshape((stateDim*nn,))))            

        if i<1:
            sol_buffer =  solve_ivp(fun=lambda t,z: jointF(t,z,par,ifOn=False), jac=lambda t,z: jointJacobian(t,z,par,ifOn=False), t_span=(times[i][0],times[i][-1]),y0=jointState, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)
        else:
            sol_buffer =  solve_ivp(fun=lambda t,z: jointF(t,z,par), jac=lambda t,z: jointJacobian(t,z,par), t_span=(times[i][0],times[i][-1]),y0=jointState, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)

        jointState = np.copy(sol_buffer.y[:,-1])

        if i<1:
            fullSol = sol_buffer.y[:,:-1]
            t = sol_buffer.t[:-1]
        else:
            fullSol = np.hstack((fullSol, sol_buffer.y[:,:]))
            t = np.hstack((t, sol_buffer.t[:]))

        trueTimeIndex = np.array([ np.sum(tSpanSimulation==t[tindex]) > 0  for tindex in range(len(t))])

        fullSol = fullSol[:,trueTimeIndex]
        t = t[trueTimeIndex]


    fullState = fullSol[:stateDim,:]
    fullSens = fullSol[stateDim:,:].reshape((stateDim, nn, len(tSpanSimulation)))

    return t, fullState, fullSens             
        
def fullSensitivities(par):    
    temp = logL(par)
    d = temp[1]
    f, sigma = par[-2:]
    _, fullState, fullSens = forwardSens(par, parIndex)
    J = f*fullState[1,ns:ne]
    x = [0]*nn     
    x[-2] = np.sum(J*(-J + d)/(f*sigma**2))
    x[-1] = np.sum((-sigma**2 + (J - d)**2)/sigma**3)
    dldcIp = np.zeros(len(J))
    for i in range(len(J)):
        val = f*(-J[i] + d[i])/sigma**2
        dldcIp[i] = val
    for ii in range(nn-2):
        x[ii] = np.sum(dldcIp * fullSens[1,ii,ns:ne])
    return x        
        
parIndex = range(len(par)) 
nn = len(parIndex)   
y = fullSensitivities(par)
print(y)