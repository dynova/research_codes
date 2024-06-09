import numpy as np
import random
import pickle
from scipy.integrate import solve_ivp
from scipy.special import loggamma, polygamma
from scipy import stats
from time import time
from numba import njit
from itertools import combinations
from math import comb
random.seed(10)

newCases = np.array([    0,     0,     0,     0,     0,     0,     0,     0,     0,
           0,     0,     0,     0,     0,     0,     0,     0,     0,
           0,     0,     0,     0,     0,     0,     0,     0,     0,
           0,     0,     0,     0,     0,     0,     0,     0,     0,
           0,     0,     0,     1,     0,     1,    10,    12,    23,
          43,    18,    41,    33,    50,   106,   108,   202,   128,
         259,   477,  1086,  1972,  3003,  3325,  5048,  6088,  5124,
        7521,  7466,  6600,  9878,  7417,  9465, 10383,  9967, 10780,
       13763, 15925, 10452, 12067, 12107, 12398, 13414, 13357, 12261,
       10765,  8549,  9942, 13544, 12026, 10046,  9250,  9076,  7474,
        7139,  7684,  9125,  9807, 11225,  8141,  4924,  4825,  6215,
        5824,  5169,  5712,  5521,  3043,  3668,  3041,  4046,  3560,
        3179,  2757,  2229,  1673,  2282,  2537,  3120,  2295,  2454,
        2263,  1612,  1838,  2216,  1840,  1387,  1942,  1624,  1268,
        1346,  2295,  1937,  1583,  1434,  1063,  1456,  1114,  1135,
        1375,  1124,   804,   789,   731,   931,   887,   826,  1042,
         688,   657,   636,   666,   677,   765,   698,   669,   590,
         610,   460,  2369,   812,   695,   722,   347,   531,   561,
         860,   918,   712,   625,   485,   563,   573,   632,   802,
         691,   708,   597,   819,   810,   704,   564,   655,   332,
         565,   689,   734,   737,   810,   789,   678,   733,   665,
         723,   649,   871,   750,   590,   549,   739,   658,   766,
         713,   763,   575,   550,   680,   804,   862,   901,   885,
         679,   539,   742,   647,   468,   626,   721,   597,   463,
         619,   586,   667,   615,   682,   645,   699,   789,   706,
         778,   895,   781,   712,   562,   571,   614,   811,   838,
         845,   653,   639,   795,   718,   967,   911,  1067,   897,
         693,   911,   752,  1014,  1026,  1259,  1179,  1116,  1097,
        1270,  1375,  1665,  1894,  1268,  1077,  1469,  1307,  2190,
        1641,  1660,  1258,  1076,  1612,  1512,  1618,  1708])

dateIndex = ['2020-01-21', '2020-01-22', '2020-01-23', '2020-01-24', '2020-01-25', '2020-01-26', '2020-01-27', '2020-01-28', '2020-01-29', '2020-01-30', '2020-01-31', '2020-02-01', '2020-02-02', '2020-02-03', '2020-02-04', '2020-02-05', '2020-02-06', '2020-02-07', '2020-02-08', '2020-02-09', '2020-02-10', '2020-02-11', '2020-02-12', '2020-02-13', '2020-02-14', '2020-02-15', '2020-02-16', '2020-02-17', '2020-02-18', '2020-02-19', '2020-02-20', '2020-02-21', '2020-02-22', '2020-02-23', '2020-02-24', '2020-02-25', '2020-02-26', '2020-02-27', '2020-02-28', '2020-02-29', '2020-03-01', '2020-03-02', '2020-03-03', '2020-03-04', '2020-03-05', '2020-03-06', '2020-03-07', '2020-03-08', '2020-03-09', '2020-03-10', '2020-03-11', '2020-03-12', '2020-03-13', '2020-03-14', '2020-03-15', '2020-03-16', '2020-03-17', '2020-03-18', '2020-03-19', '2020-03-20', '2020-03-21', '2020-03-22', '2020-03-23', '2020-03-24', '2020-03-25', '2020-03-26', '2020-03-27', '2020-03-28', '2020-03-29', '2020-03-30', '2020-03-31', '2020-04-01', '2020-04-02', '2020-04-03', '2020-04-04', '2020-04-05', '2020-04-06', '2020-04-07', '2020-04-08', '2020-04-09', '2020-04-10', '2020-04-11', '2020-04-12', '2020-04-13', '2020-04-14', '2020-04-15', '2020-04-16', '2020-04-17', '2020-04-18', '2020-04-19', '2020-04-20', '2020-04-21', '2020-04-22', '2020-04-23', '2020-04-24', '2020-04-25', '2020-04-26', '2020-04-27', '2020-04-28', '2020-04-29', '2020-04-30', '2020-05-01', '2020-05-02', '2020-05-03', '2020-05-04', '2020-05-05', '2020-05-06', '2020-05-07', '2020-05-08', '2020-05-09', '2020-05-10', '2020-05-11', '2020-05-12', '2020-05-13', '2020-05-14', '2020-05-15', '2020-05-16', '2020-05-17', '2020-05-18', '2020-05-19', '2020-05-20', '2020-05-21', '2020-05-22', '2020-05-23', '2020-05-24', '2020-05-25', '2020-05-26', '2020-05-27', '2020-05-28', '2020-05-29', '2020-05-30', '2020-05-31', '2020-06-01', '2020-06-02', '2020-06-03', '2020-06-04', '2020-06-05', '2020-06-06', '2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', '2020-06-12', '2020-06-13', '2020-06-14', '2020-06-15', '2020-06-16', '2020-06-17', '2020-06-18', '2020-06-19', '2020-06-20', '2020-06-21', '2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25', '2020-06-26', '2020-06-27', '2020-06-28', '2020-06-29', '2020-06-30', '2020-07-01', '2020-07-02', '2020-07-03', '2020-07-04', '2020-07-05', '2020-07-06', '2020-07-07', '2020-07-08', '2020-07-09', '2020-07-10', '2020-07-11', '2020-07-12', '2020-07-13', '2020-07-14', '2020-07-15', '2020-07-16', '2020-07-17', '2020-07-18', '2020-07-19', '2020-07-20', '2020-07-21', '2020-07-22', '2020-07-23', '2020-07-24', '2020-07-25', '2020-07-26', '2020-07-27', '2020-07-28', '2020-07-29', '2020-07-30', '2020-07-31', '2020-08-01', '2020-08-02', '2020-08-03', '2020-08-04', '2020-08-05', '2020-08-06', '2020-08-07', '2020-08-08', '2020-08-09', '2020-08-10', '2020-08-11', '2020-08-12', '2020-08-13', '2020-08-14', '2020-08-15', '2020-08-16', '2020-08-17', '2020-08-18', '2020-08-19', '2020-08-20', '2020-08-21', '2020-08-22', '2020-08-23', '2020-08-24', '2020-08-25', '2020-08-26', '2020-08-27', '2020-08-28', '2020-08-29', '2020-08-30', '2020-08-31', '2020-09-01', '2020-09-02', '2020-09-03', '2020-09-04', '2020-09-05', '2020-09-06', '2020-09-07', '2020-09-08', '2020-09-09', '2020-09-10', '2020-09-11', '2020-09-12', '2020-09-13', '2020-09-14', '2020-09-15', '2020-09-16', '2020-09-17', '2020-09-18', '2020-09-19', '2020-09-20', '2020-09-21', '2020-09-22', '2020-09-23', '2020-09-24', '2020-09-25', '2020-09-26', '2020-09-27', '2020-09-28', '2020-09-29', '2020-09-30', '2020-10-01', '2020-10-02', '2020-10-03', '2020-10-04', '2020-10-05', '2020-10-06', '2020-10-07', '2020-10-08', '2020-10-09', '2020-10-10', '2020-10-11', '2020-10-12', '2020-10-13', '2020-10-14', '2020-10-15',]
totalN = len(dateIndex)
tSpan = np.linspace(0, totalN-1, totalN)
tSpanSimulation = np.linspace(0, totalN, totalN+1)
t0 = 3.29475738e+01
tdelta = 2.76767105e-01
tdelta2 = 1.14291478e+02
tdelta3 = 7.90961077e+01
S0 = totalPopulation = 19216182
b     = 1.92072320e+00
mb = 0.1
relE = 1.1
relA  = 0.9
lamb  = 9.33635294e-02
fP    = 8.72255575e-01
lamb2  = 7.71149080e+00
fP2    = 7.23070977e-01
lamb3  = 6.28732210e+00
fP3    = 6.81897679e-01
kL    = 0.94
cI    = 0.12
fA    = 0.44
fH    = 0.054
kQ    = 0.0038
jQ    = 0.4
cA    = 0.26
cH    = 0.17
fR   = 0.21
fD   = 1.23472637e-01
r     = 1.52068522e+01
par = np.array([t0, tdelta, tdelta2, tdelta3, b, lamb, fP, lamb2, fP2, lamb3, fP3, S0, mb, relE, relA, kL, cI, fA, fH, kQ, jQ, cA, cH, fR, fD, r])

SM  = S0
SP  = 0
E1M = 0
E1P = 0
E2M = 0
E2P = 0
E2Q = 0
E3M = 0
E3P = 0
E3Q = 0
E4M = 0
E4P = 0
E4Q = 0
E5M = 0
E5P = 0
E5Q = 0
AM  = 0
AP  = 0
AQ  = 0
IM  = 1.
IP  = 0
IQ  = 0
IH  = 0
RM  = 0
RP  = 0
D   = 0
ccI  = 1.

IC = state = np.array([SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,ccI])
stateDim, parDim, d = len(IC), len(par), newCases

@njit
def F(t,state,par_const,ifOn=True):

	if ifOn:

		SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,ccI=state
		t0, tdelta, tdelta2, tdelta3, b, lamb, fP, S0, mb, relE, relA, kL, cI, fA, fH, kQ, jQ, cA, cH, fR, fD, r = par_const

		kfSD = lamb*fP
		krSD = lamb*(1-fP)

		return np.array([-SM*kfSD + SP*krSD - SM*b*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, SM*kfSD - SP*krSD - SP*b*mb*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, -E1M*kL - E1M*kfSD + E1P*krSD + SM*b*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, E1M*kfSD - E1P*kL - E1P*krSD + SP*b*mb*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, -E2M*kQ - E2M*kfSD + E2P*krSD + kL*(E1M - E2M), E2M*kfSD - E2P*kQ - E2P*krSD + kL*(E1P - E2P), -E2Q*kL + kQ*(E2M + E2P), -E3M*kQ - E3M*kfSD + E3P*krSD + kL*(E2M - E3M), E3M*kfSD - E3P*kQ - E3P*krSD + kL*(E2P - E3P), E2Q*kL - E3Q*kL + kQ*(E3M + E3P), -E4M*kQ - E4M*kfSD + E4P*krSD + kL*(E3M - E4M), E4M*kfSD - E4P*kQ - E4P*krSD + kL*(E3P - E4P), E3Q*kL - E4Q*kL + kQ*(E4M + E4P), -E5M*kQ - E5M*kfSD + E5P*krSD + kL*(E4M - E5M), E5M*kfSD - E5P*kQ - E5P*krSD + kL*(E4P - E5P), E4Q*kL - E5Q*kL + kQ*(E5M + E5P), -AM*cA - AM*kQ - AM*kfSD + AP*krSD + E5M*fA*kL, AM*kfSD - AP*cA - AP*kQ - AP*krSD + E5P*fA*kL, -AQ*cA + E5Q*fA*kL + kQ*(AM + AP), E5M*kL*(1 - fA) - IM*kfSD - IM*(jQ + kQ) - IM*(fH*cI + (1-fH)*cI) + IP*krSD, E5P*kL*(1 - fA) + IM*kfSD - IP*krSD - IP*(jQ + kQ) - IP*(fH*cI + (1-fH)*cI), E5Q*kL*(1 - fA) - IQ*(fH*cI + (1-fH)*cI) + (IM + IP)*(jQ + kQ), -IH*(fR*cH + (1-fR)*cH) + fH*cI*(IM + IP + IQ), AM*cA + IH*(1-fR)*cH + IM*(1-fH)*cI - RM*kfSD + RP*krSD, RM*kfSD - RP*krSD + cA*(AP + AQ) + (1-fH)*cI*(IP + IQ), IH*fR*cH, kL*(1 - fA)*(E5M + E5P) ])

	else:

		return np.zeros(np.shape(state))

def processNewImplementation(par):

    state = np.copy(IC)
    
    durations = np.array(par[:4])
    change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
    times = []

    for i in range(len(durations)+1):

        lower = change_times[i]
        upper = change_times[i+1]

        times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper))) 

    for i in range(len(times)):
        
        if i<2:
            par_const = np.array(np.hstack([par[:5],0,0,par[11:]]))
        else:
            par_const = np.array(np.hstack([par[:5],par[5+2*(i-2)],par[6+2*(i-2)],par[11:]]))
    
        if i<1:
            sol_buffer =  solve_ivp(fun=lambda t,z: F(t,z,par_const,ifOn=False), t_span=(times[i][0],times[i][-1]),y0=state, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)
        else:
            sol_buffer =  solve_ivp(fun=lambda t,z: F(t,z,par_const), t_span=(times[i][0],times[i][-1]),y0=state, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)

        state = sol_buffer.y[:,-1]

        if i==0:
            fullSol = sol_buffer.y[:,:-1]
            t = sol_buffer.t[:-1]
        elif i<len(times)-1:
            fullSol = np.hstack((fullSol, sol_buffer.y[:,:-1]))
            t = np.hstack((t, sol_buffer.t[:-1]))
        else:
            fullSol = np.hstack((fullSol, sol_buffer.y[:,:]))
            t = np.hstack((t, sol_buffer.t[:]))

        trueTimeIndex = np.array([ np.sum(tSpanSimulation==t[tindex]) > 0  for tindex in range(len(t))])

        fullSol = fullSol[:,trueTimeIndex]
        t = t[trueTimeIndex]
        
    return t, fullSol

def likelihoodFDSensNewImplementation(par):
    
    epsilon = 6.1e-6

    tSpanSimulation = np.linspace(0, len(d), len(d)+1)
    lsensFD = np.zeros((len(parIndex)))

    for j in range(len(parIndex)):

        par_p = np.copy(par)
        par_m = np.copy(par)

        perturbedInd = parIndex[j]

        diff = epsilon

        par_p[perturbedInd] += diff
        par_m[perturbedInd] -= diff
        
        durations = np.array(par_p[:4])
        change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
        times = []
        
        for i in range(len(durations)+1):

            lower = change_times[i]
            upper = change_times[i+1]

            times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper)))
        
        _,sol_p = processNewImplementation(par_p)
        
        durations = np.array(par_m[:4])
        change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
        times = []

        for i in range(len(durations)+1):

            lower = change_times[i]
            upper = change_times[i+1]

            times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper)))
           
        _,sol_m = processNewImplementation(par_m)

        Jp = par_p[-2]*(sol_p[-1,1:]-sol_p[-1,:-1])
        Jm = par_m[-2]*(sol_m[-1,1:]-sol_m[-1,:-1])
        
        pp = par_p[-1]/(par_p[-1]+Jp)
        pm = par_m[-1]/(par_m[-1]+Jm)
        
        lp = np.zeros(len(Jp))
        lm = np.zeros(len(Jm))
        
        for i in range(len(Jp)):
            if Jp[i]>0:
                lp[i] = - ((loggamma(d[i]+par_p[-1]) -loggamma(d[i]+1)-loggamma(par_p[-1])+par_p[-1]*np.log(pp[i])+d[i]*np.log(1-pp[i])))
            else:
                lp[i] = 0
                
        for i in range(len(Jm)):
            if Jm[i]>0:
                lm[i] = - ((loggamma(d[i]+par_m[-1]) -loggamma(d[i]+1)-loggamma(par_m[-1])+par_m[-1]*np.log(pm[i])+d[i]*np.log(1-pm[i])))
            else:
                lm[i] = 0                
        
        
        lsensFD[j] = np.sum(lp-lm)/2/diff
        
    return lsensFD
    
pp = len(par)
thresh = 100
cfdt = np.zeros((pp,thresh))
cfdy = [[0]*thresh]*pp
for nn in range(1,pp+1):
    ll = list(combinations(range(pp),nn))
    kk = comb(pp,nn)
    if kk < thresh and nn < pp: 
        for jj in range(kk):
            parIndex = np.array(ll[jj])
            y = likelihoodFDSensNewImplementation(par)                     
            tic = time()
            y = likelihoodFDSensNewImplementation(par)
            cfdt[nn-1,jj] = time()-tic
            cfdy[nn-1][jj] = y
    if nn == pp: 
        parIndex = np.array(ll[0])    
        for jj in range(1):
            y = likelihoodFDSensNewImplementation(par)                     
            tic = time()
            y = likelihoodFDSensNewImplementation(par)
            cfdt[nn-1,jj] = time()-tic
            cfdy[nn-1][jj] = y
    if kk > thresh and nn < pp:
        rs = random.sample(ll,thresh)    
        for jj in range(thresh):
            parIndex = np.array(rs[jj])
            y = likelihoodFDSensNewImplementation(par)         
            tic = time()
            y = likelihoodFDSensNewImplementation(par)
            cfdt[nn-1,jj] = time()-tic
            cfdy[nn-1][jj] = y           
                                   
np.save('cfdt.npy',cfdt)

with open('cfdy', 'wb') as f:
    pickle.dump(cfdy, f)