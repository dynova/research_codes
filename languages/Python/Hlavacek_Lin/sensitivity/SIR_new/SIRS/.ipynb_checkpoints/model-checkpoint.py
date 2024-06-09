#!/usr/bin/env python
# coding: utf-8

# In[82]:


import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt
newCases = np.genfromtxt('sirs.txt')[560:,1]
totalN, vv, N = len(newCases), 5000, 100
tSpan = np.linspace(0, totalN-1+vv, totalN+vv)
M1, T1, P1, M2, T2, P2, M3, T3, P3, M4, T4, P4, M5, T5, P5, gamma, t0, t1, t2, t3 = 0.0116543685, 373.349089, 91.9259976, 0.0116543685, 373.349089, 91.9259976, 0.0116543685, 373.349089, 91.9259976, 0.0116543685, 373.349089, 91.9259976, 0.0116543685, 373.349089, 91.9259976, 0.0100609532, 387.5, 387.5, 387.5, 387.5
e1,e2,e3,e4,e5,e6 = 2,1.7,2,1.25,1.8,0.1
A1, A2, A3, A4, A5, xi = 0.0464838855004532 * e1, 0.04618037483928664 * e2, 0.045733754338165046 * e3, 0.0456556809310912 * e4, 0.04465499770470104 * e5, 0.117367551 * e6
T4 = 375
par = np.array([M1, T1, A1, P1, M2, T2, A2, P2, M3, T3, A3, P3, M4, T4, A4, P4, M5, T5, A5, P5, gamma, xi, t0, t1, t2, t3])
s0, i0, r0, s, i ,r = 1 - 1/N, 1/N, 0, s0, i0, r0
IC = state = np.array([s,i,r])
stateDim, parDim = len(IC), len(par)
def F(t,state,par):
    s, i, r = state
    M1, T1, A1, P1, M2, T2, A2, P2, M3, T3, A3, P3, M4, T4, A4, P4, M5, T5, A5, P5, gamma, xi, t0, t1, t2, t3 = par  
    beta1, beta2, beta3, beta4, beta5 = M1 + A1 * np.sin(2*np.pi*t*(T1**-1) + P1), M2 + A2 * np.sin(2*np.pi*t*(T2**-1) + P2), M3 + A3 * np.sin(2*np.pi*t*(T3**-1) + P3), M4 + A4 * np.sin(2*np.pi*t*(T4**-1) + P4), M5 + A5 * np.sin(2*np.pi*t*(T5**-1) + P5)
    if t <= vv + t0:
        beta = beta1
    elif t <= vv + t0 + t1:
        beta = beta2
    elif t <= vv + t0 + t1 + t2:
        beta = beta3
    elif t <= vv + t0 + t1 + t2 + t3:
        beta = beta4
    else:
        beta = beta5
    return np.array([-beta*s*i+xi*r,beta*s*i-gamma*i,gamma*i-xi*r])
def process(par):
    sol = solve_ivp(fun=lambda t,z: F(t,z,par), t_span=(tSpan[0],tSpan[-1]),y0=state, t_eval=tSpan, method='RK45',rtol=1e-6,atol=1e-6)
    newCases = 0.08*sol.y[1,:]+0.01
    return newCases
def penalty(par):
    output = process(par)[vv:]
    Lval = np.sum((output-newCases)**2)
    return Lval        
plt.plot(tSpan[vv:]-vv,process(par)[vv:],label='model')           
plt.scatter(tSpan[vv:]-vv,newCases,s=10,color='r',label='data')
plt.legend(loc='upper right')


# In[83]:


sol = solve_ivp(fun=lambda t,z: F(t,z,par), t_span=(tSpan[0],tSpan[-1]),y0=state, t_eval=tSpan, method='RK45',rtol=1e-6,atol=1e-6)
plt.plot(tSpan[vv:]-vv,sol.y[0,vv:],label='S')           
plt.plot(tSpan[vv:]-vv,sol.y[1,vv:],label='I')           
plt.plot(tSpan[vv:]-vv,sol.y[2,vv:],label='R')           
plt.legend(loc=((0.73,0.4)))


# In[84]:


list(par)


# In[ ]:




