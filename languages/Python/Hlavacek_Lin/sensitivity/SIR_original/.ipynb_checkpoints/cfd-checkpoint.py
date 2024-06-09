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
    s, i, r = state
    t0, beta1, gamma, f, sigma = par
    if ifOn:
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
        return np.zeros(np.shape(state),np.shape(state))           

def process(par):
    sol = solve_ivp(fun=lambda t,z: F(t,z,par), t_span=(tSpan[0],tSpan[-1]),y0=state, t_eval=tSpan, method='LSODA',rtol=1e-6,atol=1e-6)
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

def processNewImplementation(par):
    state = np.copy(IC)
    durations = np.array(par[:1])
    change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
    times = []

    for i in range(len(durations)+1):

        lower = change_times[i]
        upper = change_times[i+1]

        times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper))) 

    for i in range(len(times)):
        
        if i<1:
            sol_buffer =  solve_ivp(fun=lambda t,z: F(t,z,par,ifOn=False),jac=lambda t,z: J(t,z,par,ifOn=False), t_span=(times[i][0],times[i][-1]),y0=state, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)
        else:
            sol_buffer =  solve_ivp(fun=lambda t,z: F(t,z,par,ifOn=True),jac=lambda t,z: J(t,z,par,ifOn=True), t_span=(times[i][0],times[i][-1]),y0=state, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)

        state = sol_buffer.y[:,-1]

        if i<1:
            fullSol = sol_buffer.y[:,:-1]
            t = sol_buffer.t[:-1]
        else:
            fullSol = np.hstack((fullSol, sol_buffer.y[:,:]))
            t = np.hstack((t, sol_buffer.t[:]))

        trueTimeIndex = np.array([ np.sum(tSpanSimulation==t[tindex]) > 0  for tindex in range(len(t))])

        fullSol = fullSol[:,trueTimeIndex]
        t = t[trueTimeIndex]
        
    return t, fullSol

def likelihoodFDSensNewImplementation(par):
    
    epsilon = 6.1e-6
    temp = logL(par)
    d = temp[1]
    tSpanSimulation = np.linspace(0, len(d), len(d)+1)
    lsensFD = np.zeros((len(parIndex)))

    for j in range(len(parIndex)):

        par_p = np.copy(par)
        par_m = np.copy(par)

        perturbedInd = parIndex[j]

        diff = epsilon

        par_p[perturbedInd] += diff
        par_m[perturbedInd] -= diff
        
        durations = np.array(par_p[:1])
        change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
        times = []
        
        for i in range(len(durations)+1):

            lower = change_times[i]
            upper = change_times[i+1]

            times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper)))
        
        _,sol_p = processNewImplementation(par_p)
        
        durations = np.array(par_m[:1])
        change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
        times = []

        for i in range(len(durations)+1):

            lower = change_times[i]
            upper = change_times[i+1]

            times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper)))
           
        _,sol_m = processNewImplementation(par_m)

        Jp = par_p[-2]*sol_p[1,ns:ne]
        Jm = par_m[-2]*sol_m[1,ns:ne]
               
        lp = np.zeros(len(Jp))
        lm = np.zeros(len(Jm))
        
        for i in range(len(Jp)):
            lp[i] = -0.5*np.log(2*np.pi)-np.log(par_p[-1])-0.5*(par_p[-1]**-2)*(par_p[-2]*sol_p[1,ns+i]-d[i])**2

        for i in range(len(Jm)):
            lm[i] = -0.5*np.log(2*np.pi)-np.log(par_m[-1])-0.5*(par_m[-1]**-2)*(par_m[-2]*sol_m[1,ns+i]-d[i])**2        
        
        lsensFD[j] = np.sum(lp-lm)/2/diff
        
    return lsensFD
    
parIndex = range(len(par))    
y = likelihoodFDSensNewImplementation(par)
print(y)