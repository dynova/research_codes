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
        
def dJdx(t,state,par,ifOn=True): 

    if ifOn:

        s, i, r = state
        t0, beta1, gamma, f, sigma = par  
        stateDim = len(state)  
        dflatJdx = np.array([[0, -beta, 0], [-beta, 0, 0], [0, 0, 0], [0, beta, 0], [beta, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        dJdx = np.reshape(dflatJdx, (stateDim, stateDim, stateDim))

        return dJdx
    
    else:

        return np.zeros((stateDim,stateDim,stateDim))           
        
def jointBackward(t,jointState,par,ifOn):
    
    '''RHS of the adjoint pass; first stateDim entries are the state, stateDim:2*stateDim entries are the adjoints,
    2*stateDim:2*stateDim+parDim are the sensitivities. All of which shall be integrated backward in time. '''
    
    assert len(jointState) == 2*stateDim + nn
    
    if ifOn:

        x = jointState[:stateDim]
        adj = jointState[stateDim:2*stateDim]
        
        dx = F(t,x,par,ifOn=ifOn).reshape((stateDim,))
        dadj = -(adj.reshape([1,stateDim]).dot(J(t,x,par,ifOn=ifOn))).reshape((stateDim,))
        dsens = -(adj.reshape([1,stateDim]).dot(dFdtheta_constant(t,x,par,ifOn=ifOn))).reshape((nn,)) # negative sign reflects the exchange of opper/lower limit

        return np.hstack((dx, dadj, dsens))
    
    else:
        
        return np.zeros(np.shape(jointState))
        
def jointJacobian(t,jointState,par,ifOn):

    '''The joint Jacobian of the jointBackward flow. '''
    
    assert len(jointState) == 2*stateDim + nn
    
    if ifOn:

        jointJ = np.zeros((2*stateDim + nn, 2*stateDim + nn))
        
        x = jointState[:stateDim]
        adj = jointState[stateDim:2*stateDim]

        # The Jacobian of the state is simply the Jacobian of the F
        jointJ[:stateDim, :stateDim] = J(t,x,par,ifOn=ifOn)[:,:]

        # The Jacobian of the adjoint with respect to the adjoint is simply the -Jacobian.T of the F ()
        jointJ[stateDim:2*stateDim, stateDim:2*stateDim] = - J(t,x,par,ifOn=ifOn).T[:,:]

        # The Jacobian of the adjoint with respect to the state is - adj (row vec) \cdot dJ/dx 
        dJdxval = dJdx(t,x,par,ifOn=ifOn)
        jointJ[stateDim:2*stateDim, :stateDim] = - np.tensordot(adj.reshape([1,stateDim]), dJdxval, axes=(1,0))
        
        # The Jacobian of the sens with respect to the adjoint is - dFdtheta_constant.T
        jointJ[2*stateDim:, stateDim:2*stateDim] = -dFdtheta_constant(t,x,par,ifOn=ifOn).T

        # The Jacobian of the sens with respect to the state is - adj (row vec) \cdot d2Fdthetadx_constant
        d2FdthetadxVal = d2Fdthetadx_constant(t,x,par,ifOn=ifOn)
        jointJ[2*stateDim:, :stateDim] = -np.tensordot(adj.reshape([1,stateDim]), d2FdthetadxVal, axes=(1,0))

        return jointJ

    else:

        return np.zeros((np.shape(jointState), np.shape(jointState)))      
        
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

def deltaType3(state,par,jj):            
    s, i, r = state
    t0, beta1, gamma, f, sigma = par    

    if jj == -1:
        return -F(t0,state,par,ifOn=True)
    elif jj == 0:
        beta = beta1
        return np.array([beta*s*i,-beta*s*i,0])
        
def extractDataDeltas(t, fullState, d_aug, par):

    output = np.zeros((fullState.shape[0],))
    f, sigma = par[-2:]
    
    if t >= ns and t < ne:
        alpha = fullState[1,t]
        delta = f*(-alpha*f + d_aug[t-ns])/sigma**2
    else:
        delta = 0
    
    output[1] = delta
    
    return output
    
def forwardSimulation(par, ifTease=True):

    '''Forward simulator. When ifTease=True, it returns only the snapshots of the system at the specified times,
       i.e., tSpanSimulation. When ifTease=False, it returns all the snapshots including the change times of the 
       piecewise-stationary flow. '''

    durations = np.array(par[:1])
    change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
    times = []

    for i in range(len(durations)+1):

        lower = change_times[i]
        upper = change_times[i+1]

        times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper)))

    state = np.copy(IC)

    fullSolution = []

    for i in range(len(times)):
    
        if i<1:
            sol_buffer =  solve_ivp(fun=lambda t,z: F(t,z,par,ifOn=False),jac=lambda t,z: J(t,z,par,ifOn=False), t_span=(times[i][0],times[i][-1]),y0=state, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)
        else:
            sol_buffer =  solve_ivp(fun=lambda t,z: F(t,z,par,ifOn=True),jac=lambda t,z:J(t,z,par,ifOn=True), t_span=(times[i][0],times[i][-1]),y0=state, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)

        fullSolution.append(sol_buffer.y)

        state = sol_buffer.y[:,-1]

        if i==0:
            fullState = sol_buffer.y[:,:-1]
            t = sol_buffer.t[:-1]
        else:
            fullState = np.hstack((fullState, sol_buffer.y[:,:]))
            t = np.hstack((t, sol_buffer.t[:]))

    trueTimeIndex = np.array([ np.sum(tSpanSimulation==t[tindex]) > 0  for tindex in range(len(t))])
    fullState = fullState[:,trueTimeIndex]
    t = t[trueTimeIndex]

    if ifTease:

        return t, fullState

    else:

        return times, fullSolution, t, fullState             

def adjointSensitivities(par):

    # forward pass
    times,fullSolution,t,fullState = forwardSimulation(par, ifTease=False)
    
    temp = logL(par)
    d = temp[1]
    d_aug = d

    # joint state of (x, adjoint, sensitivities)
    jointState = np.zeros(stateDim*2+nn)
    jointState[:stateDim] = fullSolution[1][:,-1]
    jointState[stateDim:2*stateDim] = extractDataDeltas(t[-1].astype('int'), fullState, d_aug, par) # terminal deltas for lambda

    dataTicks = set(tSpanSimulation)

    dynamicsTicks = []

    for i in range(len(times)):
        dynamicsTicks.append(times[i][-1])
    dynamicsTicks = set(dynamicsTicks[:-1])

    for jj in range(-1, -3, -1):

        ifOn = jj > -2
        
        time = times[jj+2]

        for i in range(len(time)-1):

            jointState[:stateDim] = fullSolution[jj+2][:,-1-i]
            sol_buffer = solve_ivp(fun=lambda t,z: jointBackward(t,z,par, ifOn=ifOn), jac = lambda t,z: jointJacobian(t,z,par, ifOn=ifOn), t_span=(times[jj+2][[-1-i,-2-i]]),y0=jointState, t_eval=times[jj+2][::-1][i:i+2], method='LSODA',rtol=1e-6,atol=1e-6)

            jointState[:] = sol_buffer.y[:,1]

            ifDataTick = time[-2-i] in dataTicks  
            ifDynamicTick = time[-2-i] in dynamicsTicks

            if ifDynamicTick:

                if ifDataTick:
                    offset = 0.5*extractDataDeltas(time[-2-i].astype(int), fullState, d_aug, par)
                    # if coincide with the data delta shocks, adjoints should be half way through the shocks
                else:
                    offset = 0

                for ii in range(1):

                    buffer = deltaType3(jointState[:stateDim], par, jj)

                    if jj>-2+parIndex[ii]:

                        # jj=-1: start to accumulate t0
                        # jj=0: start to accumulate t1
                        # jj=1: start to accumulate t2
                        # jj=2: start to accumulate t3
                        #print(buffer.shape, jointState[stateDim:2*stateDim].shape)
                        jointState[2*stateDim+ii] += np.sum(buffer*(jointState[stateDim:2*stateDim]+offset))
                        

            if ifDataTick:

                jointState[stateDim:2*stateDim] += extractDataDeltas(time[-2-i].astype(int), fullState, d_aug, par)        

    f,sigma = par[-2:]

    J = f*fullState[1,ns:ne]
    jointState[-2] = np.sum(J*(-J + d)/(f*sigma**2))
    jointState[-1] = np.sum((-sigma**2 + (J - d)**2)/sigma**3)
    
    return jointState[2*stateDim+parIndex]              
        
parIndex = np.arange(len(par)) 
nn = len(parIndex)   
y = adjointSensitivities(par)
print(y)