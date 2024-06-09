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
        
@njit
def J(t,state,par_const,ifOn=True):

    if ifOn:

        SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,ccI=state
        t0, tdelta, tdelta2, tdelta3, b, lamb, fP, S0, mb, relE, relA, kL, cI, fA, fH, kQ, jQ, cA, cH, fR, fD, r = par_const

        kfSD = lamb*fP
        krSD = lamb*(1-fP)

        return np.array([[-kfSD - b*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, krSD, 0, 0, -SM*b*relE/S0, -SM*b*mb*relE/S0, 0, -SM*b*relE/S0, -SM*b*mb*relE/S0, 0, -SM*b*relE/S0, -SM*b*mb*relE/S0, 0, -SM*b*relE/S0, -SM*b*mb*relE/S0, 0, -SM*b*relA/S0, -SM*b*mb*relA/S0, 0, -SM*b/S0, -SM*b*mb/S0, 0, 0, 0, 0, 0, 0], [kfSD, -krSD - b*mb*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, 0, 0, -SP*b*mb*relE/S0, -SP*b*mb**2*relE/S0, 0, -SP*b*mb*relE/S0, -SP*b*mb**2*relE/S0, 0, -SP*b*mb*relE/S0, -SP*b*mb**2*relE/S0, 0, -SP*b*mb*relE/S0, -SP*b*mb**2*relE/S0, 0, -SP*b*mb*relA/S0, -SP*b*mb**2*relA/S0, 0, -SP*b*mb/S0, -SP*b*mb**2/S0, 0, 0, 0, 0, 0, 0], [b*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, 0, -kL - kfSD, krSD, SM*b*relE/S0, SM*b*mb*relE/S0, 0, SM*b*relE/S0, SM*b*mb*relE/S0, 0, SM*b*relE/S0, SM*b*mb*relE/S0, 0, SM*b*relE/S0, SM*b*mb*relE/S0, 0, SM*b*relA/S0, SM*b*mb*relA/S0, 0, SM*b/S0, SM*b*mb/S0, 0, 0, 0, 0, 0, 0], [0, b*mb*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, kfSD, -kL - krSD, SP*b*mb*relE/S0, SP*b*mb**2*relE/S0, 0, SP*b*mb*relE/S0, SP*b*mb**2*relE/S0, 0, SP*b*mb*relE/S0, SP*b*mb**2*relE/S0, 0, SP*b*mb*relE/S0, SP*b*mb**2*relE/S0, 0, SP*b*mb*relA/S0, SP*b*mb**2*relA/S0, 0, SP*b*mb/S0, SP*b*mb**2/S0, 0, 0, 0, 0, 0, 0], [0, 0, kL, 0, -kL - kQ - kfSD, krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, kL, kfSD, -kL - kQ - krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, kQ, kQ, -kL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, kL, 0, 0, -kL - kQ - kfSD, krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, kL, 0, kfSD, -kL - kQ - krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, kL, kQ, kQ, -kL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, kL, 0, 0, -kL - kQ - kfSD, krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, kL, 0, kfSD, -kL - kQ - krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, kL, kQ, kQ, -kL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kL, 0, 0, -kL - kQ - kfSD, krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kL, 0, kfSD, -kL - kQ - krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kL, kQ, kQ, -kL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, fA*kL, 0, 0, -cA - kQ - kfSD, krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, fA*kL, 0, kfSD, -cA - kQ - krSD, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, fA*kL, kQ, kQ, -cA, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kL*(1 - fA), 0, 0, 0, 0, 0, -jQ - fH*cI - kQ - (1-fH)*cI - kfSD, krSD, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kL*(1 - fA), 0, 0, 0, 0, kfSD, -jQ - fH*cI - kQ - (1-fH)*cI - krSD, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kL*(1 - fA), 0, 0, 0, jQ + kQ, jQ + kQ, -fH*cI - (1-fH)*cI, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, fH*cI, fH*cI, fH*cI, -fR*cH - (1-fR)*cH, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, cA, 0, 0, (1-fH)*cI, 0, 0, (1-fR)*cH, -kfSD, krSD, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, cA, cA, 0, (1-fH)*cI, (1-fH)*cI, 0, kfSD, -krSD, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, fR*cH, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kL*(1 - fA), kL*(1 - fA), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])


    else:

        return np.zeros((len(state),len(state)))        

@njit
def dFdtheta_constant(t,state,par_const,SDstage,stateDim,parIndex,ifOn=True):

    if ifOn:

        SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,ccI=state
        t0, tdelta, tdelta2, tdelta3, b, lamb, fP, S0, mb, relE, relA, kL, cI, fA, fH, kQ, jQ, cA, cH, fR, fD, r = par_const

        kfSD = lamb*fP
        krSD = lamb*(1-fP)

        buffer = np.array([[-SM*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, -SM*fP + SP*(1 - fP), -SM*lamb - SP*lamb, SM*b*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0**2, -SM*b*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P))/S0, -SM*b*(E2M + E3M + E4M + E5M + mb*(E2P + E3P + E4P + E5P))/S0, -SM*b*(AM + AP*mb)/S0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-SP*mb*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, SM*fP - SP*(1 - fP), SM*lamb + SP*lamb, SP*b*mb*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0**2, -SP*b*mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P))/S0 - SP*b*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, -SP*b*mb*(E2M + E3M + E4M + E5M + mb*(E2P + E3P + E4P + E5P))/S0, -SP*b*mb*(AM + AP*mb)/S0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [SM*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, -E1M*fP + E1P*(1 - fP), -E1M*lamb - E1P*lamb, -SM*b*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0**2, SM*b*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P))/S0, SM*b*(E2M + E3M + E4M + E5M + mb*(E2P + E3P + E4P + E5P))/S0, SM*b*(AM + AP*mb)/S0, -E1M, 0, 0, 0, 0, 0, 0, 0, 0], [SP*mb*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, E1M*fP - E1P*(1 - fP), E1M*lamb + E1P*lamb, -SP*b*mb*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0**2, SP*b*mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P))/S0 + SP*b*(AM*relA + IM + mb*(AP*relA + IP + relE*(E2P + E3P + E4P + E5P)) + relE*(E2M + E3M + E4M + E5M))/S0, SP*b*mb*(E2M + E3M + E4M + E5M + mb*(E2P + E3P + E4P + E5P))/S0, SP*b*mb*(AM + AP*mb)/S0, -E1P, 0, 0, 0, 0, 0, 0, 0, 0], [0, -E2M*fP + E2P*(1 - fP), -E2M*lamb - E2P*lamb, 0, 0, 0, 0, E1M - E2M, 0, 0, 0, -E2M, 0, 0, 0, 0], [0, E2M*fP - E2P*(1 - fP), E2M*lamb + E2P*lamb, 0, 0, 0, 0, E1P - E2P, 0, 0, 0, -E2P, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, -E2Q, 0, 0, 0, E2M + E2P, 0, 0, 0, 0], [0, -E3M*fP + E3P*(1 - fP), -E3M*lamb - E3P*lamb, 0, 0, 0, 0, E2M - E3M, 0, 0, 0, -E3M, 0, 0, 0, 0], [0, E3M*fP - E3P*(1 - fP), E3M*lamb + E3P*lamb, 0, 0, 0, 0, E2P - E3P, 0, 0, 0, -E3P, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E2Q - E3Q, 0, 0, 0, E3M + E3P, 0, 0, 0, 0], [0, -E4M*fP + E4P*(1 - fP), -E4M*lamb - E4P*lamb, 0, 0, 0, 0, E3M - E4M, 0, 0, 0, -E4M, 0, 0, 0, 0], [0, E4M*fP - E4P*(1 - fP), E4M*lamb + E4P*lamb, 0, 0, 0, 0, E3P - E4P, 0, 0, 0, -E4P, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E3Q - E4Q, 0, 0, 0, E4M + E4P, 0, 0, 0, 0], [0, -E5M*fP + E5P*(1 - fP), -E5M*lamb - E5P*lamb, 0, 0, 0, 0, E4M - E5M, 0, 0, 0, -E5M, 0, 0, 0, 0], [0, E5M*fP - E5P*(1 - fP), E5M*lamb + E5P*lamb, 0, 0, 0, 0, E4P - E5P, 0, 0, 0, -E5P, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E4Q - E5Q, 0, 0, 0, E5M + E5P, 0, 0, 0, 0], [0, -AM*fP + AP*(1 - fP), -AM*lamb - AP*lamb, 0, 0, 0, 0, E5M*fA, 0, E5M*kL, 0, -AM, 0, -AM, 0, 0], [0, AM*fP - AP*(1 - fP), AM*lamb + AP*lamb, 0, 0, 0, 0, E5P*fA, 0, E5P*kL, 0, -AP, 0, -AP, 0, 0], [0, 0, 0, 0, 0, 0, 0, E5Q*fA, 0, E5Q*kL, 0, AM + AP, 0, -AQ, 0, 0], [0, -IM*fP + IP*(1 - fP), -IM*lamb - IP*lamb, 0, 0, 0, 0, E5M*(1 - fA), -IM, -E5M*kL, 0, -IM, -IM, 0, 0, 0], [0, IM*fP - IP*(1 - fP), IM*lamb + IP*lamb, 0, 0, 0, 0, E5P*(1 - fA), -IP, -E5P*kL, 0, -IP, -IP, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, E5Q*(1 - fA), -IQ, -E5Q*kL, 0, IM + IP, IM + IP, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, fH*(IM + IP + IQ), 0, cI*(IM + IP + IQ), 0, 0, 0, -IH, 0], [0, -RM*fP + RP*(1 - fP), -RM*lamb - RP*lamb, 0, 0, 0, 0, 0, IM*(1 - fH), 0, -IM*cI, 0, 0, AM, IH*(1 - fR), -IH*cH], [0, RM*fP - RP*(1 - fP), RM*lamb + RP*lamb, 0, 0, 0, 0, 0, (1 - fH)*(IP + IQ), 0, -cI*(IP + IQ), 0, 0, AP + AQ, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, IH*fR, IH*cH], [0, 0, 0, 0, 0, 0, 0, (1 - fA)*(E5M + E5P), 0, -kL*(E5M + E5P), 0, 0, 0, 0, 0, 0]])        

        output = np.zeros((stateDim,nn))

        for jj in range(nn):
            if parIndex[jj] in [0,1,2,3,24,25]:
                pass
            elif parIndex[jj] == 4:
                output[:,jj] = buffer[:,0]  
            elif parIndex[jj] in [5,6,7,8,9,10]:
                if SDstage>=0:
                    if parIndex[jj] == 5+2*SDstage:
                        output[:,jj] = buffer[:,1]
                    elif parIndex[jj] == 6+2*SDstage:
                        output[:,jj] = buffer[:,2]
            for ii in range(11,24):                    
                if parIndex[jj] == ii:
                    output[:,jj] = buffer[:,ii-8]           
        return output
    else:     
        return np.zeros((stateDim,nn))                                        

def jointF(t,jointState,par_const,SDstage,stateDim,parIndex,ifOn=True):

    assert len(jointState)==stateDim+stateDim*nn

    if ifOn:

        x = jointState[:stateDim]
        s = jointState[stateDim:].reshape((stateDim,nn))

        dx = F(t,x,par_const,ifOn)
        ds = (J(t,x,par_const,ifOn).dot(s)+ dFdtheta_constant(t,x,par_const,SDstage,stateDim,parIndex,ifOn)).reshape((stateDim*nn,))

        return np.hstack((dx,ds))

    else:

        return np.zeros(np.shape(jointState))
        
def deltaType3(state,par,SDstage):            

    SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,ccI=state
    t0, tdelta, tdelta2, tdelta3, b, lamb, fP, lamb2, fP2, lamb3, fP3, S0, mb, relE, relA, kL, cI, fA, fH, kQ, jQ, cA, cH, fR, fD, r = par

    if SDstage==-2:

        raise ValueError("There is no delta associated with the onset/stage=0!")

    elif SDstage==-1:

        par_const = np.array(np.hstack([par[:5],0,0,par[11:]]))         

        return -F(t0,state,par_const,ifOn=True)

    elif SDstage==0:

        dkf = par[5]*par[6]
        dkr = par[5]*(1-par[6])

    else:

        dkf = par[5+2*SDstage]*par[6+2*SDstage]-par[5+2*(SDstage-1)]*par[6+2*(SDstage-1)]
        dkr = par[5+2*SDstage]*(1.-par[6+2*SDstage])-par[5+2*(SDstage-1)]*(1.-par[6+2*(SDstage-1)])

    return np.array([[SM*dkf - SP*dkr], [-SM*dkf + SP*dkr], [E1M*dkf - E1P*dkr], [-E1M*dkf + E1P*dkr], [E2M*dkf - E2P*dkr], [-E2M*dkf + E2P*dkr], [0], [E3M*dkf - E3P*dkr], [-E3M*dkf + E3P*dkr], [0], [E4M*dkf - E4P*dkr], [-E4M*dkf + E4P*dkr], [0], [E5M*dkf - E5P*dkr], [-E5M*dkf + E5P*dkr], [0], [AM*dkf - AP*dkr], [-AM*dkf + AP*dkr], [0], [IM*dkf - IP*dkr], [-IM*dkf + IP*dkr], [0], [0], [RM*dkf - RP*dkr], [-RM*dkf + RP*dkr], [0], [0]]).reshape((stateDim,))

def forwardSens(par, stateDim, parIndex, tSpanSimulation, IC):

    jointState = np.zeros(stateDim*(nn+1))
    jointState[:stateDim] = IC[:]

    durations = np.array(par[:4])
    change_times = np.hstack([tSpanSimulation[0], np.cumsum(durations), tSpanSimulation[-1]])
    times = []

    for i in range(len(durations)+1):

        lower = change_times[i]
        upper = change_times[i+1]

        times.append(np.hstack((lower, tSpanSimulation[ np.logical_and(tSpanSimulation > lower, tSpanSimulation < upper)], upper)))


    for i in range(len(times)):

        SDstage = i-2 # first two phases are 1. no disease, 2. no SD. 

        if SDstage<0:
            par_const = np.array(np.hstack([par[:5],0,0,par[11:]]))
        else:
            par_const = np.array(np.hstack([par[:5],par[5+2*SDstage],par[6+2*SDstage],par[11:]]))

        ## adding the delta shocks

        for ii in range(nn):    
            if SDstage>-2+parIndex[ii]:
                # SDstage=-1: start to accumulate t0
                # SDstage=0: start to accumulate t1
                # SDstage=1: start to accumulate t2
                # SDstage=2: start to accumulate t3
                x,s = jointState[:stateDim],jointState[stateDim:].reshape((stateDim, nn))
                buffer = deltaType3(x, par, SDstage)
                s[:,ii] += buffer    
                jointState = np.hstack((x, s.reshape((stateDim*nn,))))

        if i<1:
            sol_buffer =  solve_ivp(fun=lambda t,z: jointF(t,z,par_const,SDstage=SDstage,stateDim=stateDim,parIndex=parIndex,ifOn=False), t_span=(times[i][0],times[i][-1]),y0=jointState, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)
        else:
            sol_buffer =  solve_ivp(fun=lambda t,z: jointF(t,z,par_const,SDstage=SDstage,stateDim=stateDim,parIndex=parIndex), t_span=(times[i][0],times[i][-1]),y0=jointState, t_eval=times[i], method='LSODA',rtol=1e-6,atol=1e-6)

        jointState = np.copy(sol_buffer.y[:,-1])

        if i==0:
            fullSol = sol_buffer.y[:,:-1]
            t = sol_buffer.t[:-1]
        elif i<len(times)-2:
            fullSol = np.hstack((fullSol, sol_buffer.y[:,:-1]))
            t = np.hstack((t, sol_buffer.t[:-1]))
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
    d = newCases
    tSpanSimulation = np.linspace(0, len(d), len(d)+1)
    stateDim = len(IC)
    fD,r = par[-2:]
    
    if nn > 1:
        if parIndex[-2] == 24 and parIndex[-1] == 25:
            _, fullState, fullSens = forwardSens(par, stateDim, parIndex, tSpanSimulation, IC)
            J = fD*(fullState[-1,1:]-fullState[-1,:-1])
            dldcIp = np.zeros(len(J))
            for i in range(len(J)):
                if J[i]!=0:
                    val = fD*r*(J[i] - d[i])/(J[i]*(J[i] + r))
                else:
                    val = 0.
                dldcIp[i] = val
            x = [0]*nn
            x[-2] = np.sum(r*(J - d)/(fD*(J + r)))
            x[-1] = np.sum((-J + d + (J + r)*(-np.log(r/(J + r)) + polygamma(0, r) - polygamma(0, d + r)))/(J + r))                
            for ii in range(nn-2):
                x[ii] = np.sum(dldcIp * (fullSens[-1,ii,1:] - fullSens[-1,ii,:-1]))
            return x
        elif parIndex[-1] == 24:
            _, fullState, fullSens = forwardSens(par, stateDim, parIndex, tSpanSimulation, IC)
            J = fD*(fullState[-1,1:]-fullState[-1,:-1])
            dldcIp = np.zeros(len(J))
            for i in range(len(J)):
                if J[i]!=0:
                    val = fD*r*(J[i] - d[i])/(J[i]*(J[i] + r))
                else:
                    val = 0.
                dldcIp[i] = val
            x = [0]*nn
            x[-1] = np.sum(r*(J - d)/(fD*(J + r)))
            for ii in range(nn-1):
                x[ii] = np.sum(dldcIp * (fullSens[-1,ii,1:] - fullSens[-1,ii,:-1]))
            return x
        elif parIndex[-1] == 25:
            _, fullState, fullSens = forwardSens(par, stateDim, parIndex, tSpanSimulation, IC)
            J = fD*(fullState[-1,1:]-fullState[-1,:-1])
            dldcIp = np.zeros(len(J))
            for i in range(len(J)):
                if J[i]!=0:
                    val = fD*r*(J[i] - d[i])/(J[i]*(J[i] + r))
                else:
                    val = 0.
                dldcIp[i] = val
            x = [0]*nn
            x[-1] = np.sum((-J + d + (J + r)*(-np.log(r/(J + r)) + polygamma(0, r) - polygamma(0, d + r)))/(J + r))
            for ii in range(nn-1):
                x[ii] = np.sum(dldcIp * (fullSens[-1,ii,1:] - fullSens[-1,ii,:-1]))
            return x      
        else:
            _, fullState, fullSens = forwardSens(par, stateDim, parIndex, tSpanSimulation, IC)
            J = fD*(fullState[-1,1:]-fullState[-1,:-1])
            dldcIp = np.zeros(len(J))
            for i in range(len(J)):
                if J[i]!=0:
                    val = fD*r*(J[i] - d[i])/(J[i]*(J[i] + r))
                else:
                    val = 0.
                dldcIp[i] = val
            x = [0]*nn
            for ii in range(nn):
                x[ii] = np.sum(dldcIp * (fullSens[-1,ii,1:] - fullSens[-1,ii,:-1]))
            return x    
    else:
        if parIndex[-1] == 24:
            _, fullState, fullSens = forwardSens(par, stateDim, parIndex, tSpanSimulation, IC)
            J = fD*(fullState[-1,1:]-fullState[-1,:-1])
            x = [0]*nn
            x[-1] = np.sum(r*(J - d)/(fD*(J + r)))
            return x
        elif parIndex[-1] == 25:
            _, fullState, fullSens = forwardSens(par, stateDim, parIndex, tSpanSimulation, IC)
            J = fD*(fullState[-1,1:]-fullState[-1,:-1])
            x = [0]*nn
            x[-1] = np.sum((-J + d + (J + r)*(-np.log(r/(J + r)) + polygamma(0, r) - polygamma(0, d + r)))/(J + r))
            return x      
        else:
            x = [0]*nn
            _, fullState, fullSens = forwardSens(par, stateDim, parIndex, tSpanSimulation, IC)
            J = fD*(fullState[-1,1:]-fullState[-1,:-1])
            dldcIp = np.zeros(len(J))
            for i in range(len(J)):
                if J[i]!=0:
                    val = fD*r*(J[i] - d[i])/(J[i]*(J[i] + r))
                else:
                    val = 0.
                dldcIp[i] = val
            x[-1] = np.sum(dldcIp * (fullSens[-1,0,1:] - fullSens[-1,0,:-1]))
            return x        

pp = len(par)
thresh = 100
fwdt = np.zeros((pp,thresh))
fwdy = [[0]*thresh]*pp
for nn in range(1,pp+1):
    ll = list(combinations(range(pp),nn))
    kk = comb(pp,nn)
    if kk < thresh and nn < pp: 
        for jj in range(kk):
            parIndex = np.array(ll[jj])
            y = np.array(fullSensitivities(par))                     
            tic = time()
            y = np.array(fullSensitivities(par))
            fwdt[nn-1,jj] = time()-tic
            fwdy[nn-1][jj] = y
    if nn == pp: 
        parIndex = np.array(ll[0])    
        for jj in range(1):
            y = np.array(fullSensitivities(par))                     
            tic = time()
            y = np.array(fullSensitivities(par))
            fwdt[nn-1,jj] = time()-tic
            fwdy[nn-1][jj] = y
    if kk > thresh and nn < pp:
        rs = random.sample(ll,thresh)    
        for jj in range(thresh):
            parIndex = np.array(rs[jj])
            y = np.array(fullSensitivities(par))         
            tic = time()
            y = np.array(fullSensitivities(par))
            fwdt[nn-1,jj] = time()-tic
            fwdy[nn-1][jj] = y             
                                   
np.save('fwdt.npy',fwdt)

with open('fwdy', 'wb') as f:
    pickle.dump(fwdy, f)