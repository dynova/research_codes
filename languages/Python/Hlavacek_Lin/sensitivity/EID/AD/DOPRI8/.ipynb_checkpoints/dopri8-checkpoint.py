import numpy as np
import diffrax
import jax
import jax.numpy as jnp
import jax.scipy.stats as stats
from jax.scipy.special import gammaln
from time import time
from jax import grad
jax.config.update("jax_enable_x64", True)

newCases = jnp.array([    0,     0,     0,     0,     0,     0,     0,     0,     0,
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
tSpan = jnp.linspace(0, totalN, totalN+1)
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
par = jnp.array([t0, tdelta, tdelta2, tdelta3, b, lamb, fP, lamb2, fP2, lamb3, fP3, S0, mb, relE, relA, kL, cI, fA, fH, kQ, jQ, cA, cH, fR, fD, r])

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

IC = state = jnp.array([SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,ccI])

@jax.jit
def RHS(t,state,par):
    SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,ccI = state
    t0, tdelta, tdelta2, tdelta3, b, lamb, fP, lamb2, fP2, lamb3, fP3, S0, mb, relE, relA, kL, cI, fA, fH, kQ, jQ, cA, cH, fR, fD, r = par
    kfSD = (lamb*fP)*jnp.heaviside(t-t0-tdelta,0.5)+(lamb2*fP2-lamb*fP)*jnp.heaviside(t-t0-tdelta-tdelta2,0.5)+(lamb3*fP3-lamb2*fP2)*jnp.heaviside(t-t0-tdelta-tdelta2-tdelta3,0.5)
    krSD = (lamb*(1-fP))*jnp.heaviside(t-t0-tdelta,0.5)+(lamb2*(1-fP2)-lamb*(1-fP))*jnp.heaviside(t-t0-tdelta-tdelta2,0.5)+(lamb3*(1-fP3)-lamb2*(1-fP2))*jnp.heaviside(t-t0-tdelta-tdelta2-tdelta3,0.5)
    totalInfectious = IM + relE*(E2M+E3M+E4M+E5M)+relA*AM + mb*(IP + relE*(E2P+E3P+E4P+E5P)+relA*AP)
    fSM = b*SM*totalInfectious/S0
    fSP = mb*b*SP*totalInfectious/S0
    return jnp.heaviside(t-t0,0.5)*jnp.array([-fSM - (kfSD*SM-krSD*SP), -fSP + (kfSD*SM-krSD*SP), fSM - kL*E1M-(kfSD*E1M-krSD*E1P), fSP - kL*E1P+(kfSD*E1M-krSD*E1P),kL*(E1M-E2M)-kQ*E2M-(kfSD*E2M-krSD*E2P),kL*(E1P-E2P)-kQ*E2P+(kfSD*E2M-krSD*E2P),kQ*(E2M+E2P)-kL*E2Q,kL*(E2M-E3M)-kQ*E3M-(kfSD*E3M-krSD*E3P),kL*(E2P-E3P)-kQ*E3P+(kfSD*E3M-krSD*E3P),kQ*(E3M+E3P)-kL*E3Q + kL*E2Q,kL*(E3M-E4M)-kQ*E4M-(kfSD*E4M-krSD*E4P),kL*(E3P-E4P)-kQ*E4P+(kfSD*E4M-krSD*E4P),kQ*(E4M+E4P)-kL*E4Q + kL*E3Q,kL*(E4M-E5M)-kQ*E5M-(kfSD*E5M-krSD*E5P),kL*(E4P-E5P)-kQ*E5P+(kfSD*E5M-krSD*E5P),kQ*(E5M+E5P)-kL*E5Q + kL*E4Q,fA*kL*E5M-kQ*AM-(kfSD*AM-krSD*AP)-cA*AM,fA*kL*E5P-kQ*AP+(kfSD*AM-krSD*AP)-cA*AP,fA*kL*E5Q+kQ*(AM+AP)-cA*AQ,(1-fA)*kL*E5M-(kQ+jQ)*IM-(kfSD*IM-krSD*IP)-((fH*cI)+((1-fH)*cI))*IM,(1-fA)*kL*E5P-(kQ+jQ)*IP+(kfSD*IM-krSD*IP)-((fH*cI)+((1-fH)*cI))*IP,(1-fA)*kL*E5Q+(kQ+jQ)*(IP+IM)-((fH*cI)+((1-fH)*cI))*IQ,(fH*cI)*(IM+IP+IQ)-(((1-fR)*cH)+(fR*cH))*IH,cA*AM+((1-fH)*cI)*IM+((1-fR)*cH)*IH-(kfSD*RM-krSD*RP),cA*(AP+AQ)+((1-fH)*cI)*(IP+IQ)+(kfSD*RM-krSD*RP),(fR*cH)*IH,(1-fA)*kL*(E5M+E5P)])

@jax.jit
def logL(par):
    terms = diffrax.ODETerm(RHS)
    solver = diffrax.Dopri8()
    t0 = tSpan[0]
    t1 = tSpan[-1]
    dt0 = 1e-4
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
    args = (par),
    saveat=saveat,
    stepsize_controller=stepsize_controller,
    max_steps = int(1e8),
    )
    output = par[-2]*(sol.ys[1:,-1]-sol.ys[:-1,-1])
    logLval = 0
    
    for i in range(len(newCases)):
        r = par[-1]
        prob = jnp.clip(r/(r+output[i]), 1e-10, 1-1e-10)
        logLval = logLval + gammaln(newCases[i]+r)-gammaln(newCases[i]+1)-gammaln(r)+r*jnp.log(prob)+newCases[i]*jnp.log(1-prob)
    return -logLval

count = -1
tt = np.zeros(7)
ss = np.zeros((len(par),len(tt)))
for tolerance in jnp.logspace(-8,-2,len(tt)):
    count = count + 1
    tic = time()
    ss[:,count] = grad(logL)(par)
    tt[count] = time()-tic

np.save('tt.npy',tt)
np.save('ss.npy',ss)