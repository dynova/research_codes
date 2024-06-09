from numpy import *
from scipy.stats import loguniform
from scipy.integrate import solve_ivp
from scipy.special import loggamma
from scipy.optimize import minimize
from sys import argv
from numba import njit
from matplotlib import cm
import random
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':24})
colors = cm.plasma(linspace(0,1,12))
random.seed(10)

nStages = 4
nMutations = 2
locationIndex = 0
nms = ['NYC','LA','Chicago','Dallas','Houston','DC','Philadelphia','Miami','Atlanta','Boston','Phoenix','SanFrancisco','Riverside','Detroit','Seattle']
population = [19216182,13214799,9458539,7573136,7066141,6280487,6166488,6102434,6020364,4948203,4873019,4731803,4650631,4319629,3979845]
prefix = nms[locationIndex]
totalPopulation = population[locationIndex]
filePrefix = prefix+'-n'+str(nStages)
NYTNewCases = array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,10,12,23,43,18,41,33,50,106,108,202,128,259,477,1086,1972,3003,3325,5048,6088,5124,7521,7466,6600,9878,7417,9465,10383,9967,10780,13763,15925,10452,12067,12107,12398,13414,13357,12261,10765,8549,9942,13544,12026,10046,9250,9076,7474,7139,7684,9125,9807,11225,8141,4924,4825,6215,5824,5169,5712,5521,3043,3668,3041,4046,3560,3179,2757,2229,1673,2282,2537,3120,2295,2454,2263,1612,1838,2216,1840,1387,1942,1624,1268,1346,2295,1937,1583,1434,1063,1456,1114,1135,1375,1124,804,789,731,931,887,826,1042,688,657,636,666,677,765,698,669,590,610,460,2369,812,695,722,347,531,561,860,918,712,625,485,563,573,632,802,691,708,597,819,810,704,564,655,332,565,689,734,737,810,789,678,733,665,723,649,871,750,590,549,739,658,766,713,763,575,550,680,804,862,901,885,679,539,742,647,468,626,721,597,463,619,586,667,615,682,645,699,789,706,778,895,781,712,562,571,614,811,838,845,653,639,795,718,967,911,1067,897,693,911,752,1014,1026,1259,1179,1116,1097,1270,1375,1665,1894,1268,1077,1469,1307,2190,1641,1660,1258,1076,1612,1512,1618,1708,1969,1891,1573,1479,2187,1774,1770,2733,1933,1637,2561,2673,2883,3040,2193,2607,2086,2968,3091,3208,3547,4440,3728,3452,5143,5258,5459,5623,6578,5712,3798,6313,6199,6196,5944,6724,6150,6031,6186,6384,7666,7761,5351,8066,6863,8009,9123,9600,10827,10043,10323,8457,10067,9899,9947,8685,11351,9204,9181,9298,10060,9167,10515,8785,9165,8537,9931,11006,12090,11826,10083,7520,8300,10852,11582,13406,14790,14473,10502,42774,13515,13996,17302,17032,17188,14730,13976,13162,16083,15118,17614,15227,13397,12712,12369,11756,13021,14002,14894,12617,12523,11580,10920,13926,14070,13568,11560,9789,8723,7026,6742,8452,12313,11241,9144,9836,8004,10679,9826,9941,8000,6338,8358,7879,10449,9397,8639,6154,7206,8153,7348,9656,9313,9571,7888,7921,7970,8021,8464,9334,9559,8356,6194,9157,7868,8802,10123,7970,7998,8296,8547,8549,7373,7066,8754,4540,6291,8953,21879,9749,10767,10385,10632,10029,7924,8799,11309,9671,9111,8718,6933,7877,7904,9998,9384,8018,7199,6063,6440,6712,7621,8022,6979,6131,5530,5520,5543,5904,5260,4843,3904,-5948,3934,3333,3852,3663,3520,2594,2220,1984,3625,2570,2026,2028,1447,1101,1152,1560,1421,1219,1221,1374,996,1537,999,1314,1400,962,781,536,849,697,989,758,644,661,448,448,425,514,652,768,439,481,498,435,675,739,494,411,459,462,412,499,547,444,337,333,380,383,520,531,446,407,363,361,326,474,602,352,337,388,363,478,823,882,852,944,722,835,854,1104,831,1695,1233,1044,1422,1390,1963,2086,2303,2088,1753,1862,2098,2884,2879,3363,2649,2489,3331,2997,3688,3890,4257,3534,1468,1910,7266,4872,4446,5295,3662,3979,3822,4300,5033,4395,6052,4087,3597,3822,3812,5163,5360,5026,2296,5328,3945,4075,4889,4714,2648,6101,3288,2102,3342,5757,5362,5327,2443,4739,4267,4430,5426,5223,2831,6601,3398,4550,4230,4446,4408,2638,5462,1460,4692,3526,4119,3864,2261,1783,5768,3141,4342,4293,3675,2193,1497,5579,2858,3102,3062,3982,1970,1381,4025,2754,2698,2829,2860,1368,1106,3827,2309,2172,2148,3032,1500,1423])
vr_daily=load('/Users/amallela/Documents/bmab/vac.npy')
x1 = linspace(1, len(vr_daily[locationIndex])+2, len(vr_daily[locationIndex])+2) 
y1 = vr_daily[locationIndex]
y1 = append(y1, [0,0])

@njit
def fv(x1,y1,t):
    xm=mean(x1)
    ym=mean(y1)
    sumnr=0
    sumdr=0
    length=len(x1)
    for i in range(0,length):
        sumnr=sumnr+((x1[i]-xm)*(y1[i]-ym))
        sumdr=sumdr+((x1[i]-xm)*(x1[i]-xm))
    m=sumnr/sumdr
    c=ym-(m*xm)
    return((m*t)+c)
    
tSpan = linspace(0, len(NYTNewCases)-1, len(NYTNewCases))    
tSpanSimulation = linspace(0, len(NYTNewCases), len(NYTNewCases)+1)    

S0 = totalPopulation
mb = 0.1
relEI = 1.1
relA = 0.9
kL = 0.94
gI = 0.12
fA = 0.44
fH = 0.054
kH = fH*gI
kR = (1-fH)*gI
kQ = 0.0038
jQ = 0.4
gA = 0.26
gH = 0.17
CFR = 0.21
kHD = CFR*gH
kHR = (1-CFR)*gH
kV = 0.3
kR25 = (1-fH/25)*gI
f0= 0.9
f1= 0.81
f2= 0.69
     
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
cI  = 1.
V1  = 0
V2  = 0
V3  = 0
V4  = 0
V5  = 0
V6  = 0
RV  = 0
SV1 = 0
SV2 = 0
SV3 = 0
SV4 = 0
EV = 0
AV = 0
IV = 0
HV = 0

state = [SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,V1,V2,V3,V4,V5,V6,SV1,SV2,SV3,SV4,EV,AV,IV,HV,RV,cI]
state = array(state)

@njit
def RHS(t,state,par):

    SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,V1,V2,V3,V4,V5,V6,SV1,SV2,SV3,SV4,EV,AV,IV,HV,RV,cI= state
    
    ts = par[:nStages+1]
    lambs = par[nStages+1:2*nStages+1]
    fPs = par[2*nStages+1:3*nStages+1]
    tmutation = par[3*nStages+1:3*nStages+1+nMutations]
    mutationMultiplier = par[3*nStages+1+nMutations:3*nStages+1+2*nMutations]
    
    b, fDI, r = par[nStages*3+2*nMutations+1:nStages*3+2*nMutations+4]
    
    cumsumts = cumsum(ts)[1:]
    cumsumts = append(cumsumts,inf)
    stage = where(cumsumts>t)[0][0]-1
    
    cumsumtmutation = cumsum(tmutation)[:]
    cumsumtmutation = append(cumsumtmutation,inf)
    mutationIndex = where(cumsumtmutation>t)[0][0]-1
   
    kfSD = (stage==-1)*0+(stage>=0)*lambs[stage]*fPs[stage]
    krSD = (stage==-1)*0+(stage>=0)*lambs[stage]*(1.-fPs[stage])
    
    b = (mutationIndex==-1)*b + (mutationIndex>=0)*b*mutationMultiplier[mutationIndex]
    hmutation1 = (cumsumtmutation[0]<t)*1.0
    hmutation2 = (cumsumtmutation[1]<t)*1.0
    
    totalInfectious = IM + relEI*(E2M+E3M+E4M+E5M+EV)+relA*(AM+AV) + mb*(IP + relEI*(E2P+E3P+E4P+E5P)+relA*AP)+IV
    fSM = b*SM*totalInfectious/S0
    fSP = mb*b*SP*totalInfectious/S0
    totalVacEligible=SM+SP+RM+RP+E1P+E2P+E3P+E4P+E5P+E1M+E2M+E3M+E4M+E5M+AM+AP
    mu = fv(x1,y1,t)

    return (t >= ts[0])*array([-fSM - (kfSD*SM-krSD*SP) - mu*S0*SM/totalVacEligible,    #SM
                  -fSP + (kfSD*SM-krSD*SP) - mu*S0*SP/totalVacEligible,    #SP
                  fSM - kL*E1M-(kfSD*E1M-krSD*E1P), #E1M
                  fSP - kL*E1P+(kfSD*E1M-krSD*E1P),     #E1P
                  kL*(E1M-E2M)-kQ*E2M-(kfSD*E2M-krSD*E2P),  #E2M
                  kL*(E1P-E2P)-kQ*E2P+(kfSD*E2M-krSD*E2P),  #E2P
                  kQ*(E2M+E2P)-kL*E2Q,                      #E2Q
                  kL*(E2M-E3M)-kQ*E3M-(kfSD*E3M-krSD*E3P),  #E3M
                  kL*(E2P-E3P)-kQ*E3P+(kfSD*E3M-krSD*E3P),  #E3P
                  kQ*(E3M+E3P)-kL*E3Q + kL*E2Q,             #E3Q
                  kL*(E3M-E4M)-kQ*E4M-(kfSD*E4M-krSD*E4P),  #E4M
                  kL*(E3P-E4P)-kQ*E4P+(kfSD*E4M-krSD*E4P),  #E4P
                  kQ*(E4M+E4P)-kL*E4Q + kL*E3Q,             #E4Q
                  kL*(E4M-E5M)-kQ*E5M-(kfSD*E5M-krSD*E5P),  #E5M
                  kL*(E4P-E5P)-kQ*E5P+(kfSD*E5M-krSD*E5P),  #E5P
                  kQ*(E5M+E5P)-kL*E5Q + kL*E4Q,             #E5Q
                  fA*kL*E5M-kQ*AM-(kfSD*AM-krSD*AP)-gA*AM,  #AM
                  fA*kL*E5P-kQ*AP+(kfSD*AM-krSD*AP)-gA*AP,  #AP
                  fA*kL*E5Q+kQ*(AM+AP)-gA*AQ,               #AQ
                  (1-fA)*kL*E5M-(kQ+jQ)*IM-(kfSD*IM-krSD*IP)-(kH+kR)*IM,    #IM
                  (1-fA)*kL*E5P-(kQ+jQ)*IP+(kfSD*IM-krSD*IP)-(kH+kR)*IP,    #IP
                  (1-fA)*kL*E5Q+(kQ+jQ)*(IP+IM)-(kH+kR)*IQ, #IQ
                  kH*(IM+IP+IQ)-(kHR+kHD)*IH,   #IH
                  gA*AM+kR*IM+kHR*IH-(kfSD*RM-krSD*RP)- mu*S0*RM/totalVacEligible,  #RM
                  gA*(AP+AQ)+kR*(IP+IQ)+(kfSD*RM-krSD*RP)- mu*S0*RP/totalVacEligible,#RP
                  kHD*(IH+HV), #D
                  mu*S0*(SM+SP)/totalVacEligible-kV*V1-b*V1*totalInfectious/S0,  #V1
                  kV*(V1-V2)-b*V2*totalInfectious/S0, #V2
                  kV*(V2-V3)-b*V3*totalInfectious/S0, #V3
                  kV*(V3-V4)-b*V4*totalInfectious/S0, #V4
                  kV*(V4-V5)-b*V5*totalInfectious/S0, #V5
                  kV*(V5-V6)-b*V6*totalInfectious/S0, #V6
                  kV*(1-f0)*V6-b*SV1*totalInfectious/S0, #SV1
                  kV*(f0-f1)*V6-b*hmutation1*SV2*totalInfectious/S0, #SV2
                  kV*(f1-f2)*V6-b*hmutation2*SV3*totalInfectious/S0, #SV3
                  kV*f2*V6, #SV4
                  b*(V1+V2+V3+V4+V5+V6+SV1+hmutation1*SV2+hmutation2*SV3)*totalInfectious/S0-kL*EV/5, #EV
                  fA*kL*EV/5-gA*AV, #AV
                  (1-fA)*kL*EV/5-(kH+kR)*IV, #IV
                  kH*IV/25-(kHR+kHD)*HV, #HV
                  mu*S0*(RM+RP)/totalVacEligible+gA*AV+kR25*IV+kHR*HV, #RV 
                  (1-fA)*kL*(E5M+E5P+E5Q)+(1-fA)*kL*EV/5])  #cI              

@njit
def Jacobian(t,state,par):

    SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,V1,V2,V3,V4,V5,V6,SV1,SV2,SV3,SV4,EV,AV,IV,HV,RV,cI= state
    
    ts = par[:nStages+1]
    lambs = par[nStages+1:2*nStages+1]
    fPs = par[2*nStages+1:3*nStages+1]
    tmutation = par[3*nStages+1:3*nStages+1+nMutations]
    mutationMultiplier = par[3*nStages+1+nMutations:3*nStages+1+2*nMutations]
    
    b, fDI, r = par[nStages*3+2*nMutations+1:nStages*3+2*nMutations+4]
    
    cumsumts = cumsum(ts)[1:]
    cumsumts = append(cumsumts,inf)
    stage = where(cumsumts>t)[0][0]-1
    
    cumsumtmutation = cumsum(tmutation)[:]
    cumsumtmutation = append(cumsumtmutation,inf)
    mutationIndex = where(cumsumtmutation>t)[0][0]-1
   
    kfSD = (stage==-1)*0+(stage>=0)*lambs[stage]*fPs[stage]
    krSD = (stage==-1)*0+(stage>=0)*lambs[stage]*(1.-fPs[stage])
    
    b = (mutationIndex==-1)*b + (mutationIndex>=0)*b*mutationMultiplier[mutationIndex]
    hmutation1 = (cumsumtmutation[0]<t)*1.0 #Heaviside function
    hmutation2 = (cumsumtmutation[1]<t)*1.0 #Heaviside function
    
    totalInfectious = IM + relEI*(E2M+E3M+E4M+E5M+EV)+relA*(AM+AV) + mb*(IP + relEI*(E2P+E3P+E4P+E5P)+relA*AP)+IV
    fSM = b*SM*totalInfectious/S0
    fSP = mb*b*SP*totalInfectious/S0
    totalVacEligible=SM+SP+RM+RP+E1P+E2P+E3P+E4P+E5P+E1M+E2M+E3M+E4M+E5M+AM+AP
    mu = fv(x1,y1,t)

    return (t >= ts[0])*array([[S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP) - kfSD - b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + krSD,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*relEI/S0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*mb*relEI/S0,0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*relEI/S0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*mb*relEI/S0,0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*relEI/S0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*mb*relEI/S0,0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*relEI/S0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*mb*relEI/S0,0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*relA/S0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - SM*b*mb*relA/S0,0,-SM*b/S0,-SM*b*mb/S0,0,0,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,S0*SM*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,0,0,0,0,0,0,0,0,0,0,-SM*b*relEI/S0,-SM*b*relA/S0,-SM*b/S0,0,0,0,],[S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + kfSD,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP) - krSD + b*mb*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb*relEI/S0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb**2*relEI/S0,0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb*relEI/S0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb**2*relEI/S0,0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb*relEI/S0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb**2*relEI/S0,0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb*relEI/S0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb**2*relEI/S0,0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb*relA/S0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + SP*b*mb**2*relA/S0,0,SP*b*mb/S0,SP*b*mb**2/S0,0,0,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,S0*SP*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,0,0,0,0,0,0,0,0,0,0,SP*b*mb*relEI/S0,SP*b*mb*relA/S0,SP*b*mb/S0,0,0,0,],[b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,-kL - kfSD,krSD,SM*b*relEI/S0,SM*b*mb*relEI/S0,0,SM*b*relEI/S0,SM*b*mb*relEI/S0,0,SM*b*relEI/S0,SM*b*mb*relEI/S0,0,SM*b*relEI/S0,SM*b*mb*relEI/S0,0,SM*b*relA/S0,SM*b*mb*relA/S0,0,SM*b/S0,SM*b*mb/S0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,SM*b*relEI/S0,SM*b*relA/S0,SM*b/S0,0,0,0,],[0,b*mb*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,kfSD,-kL - krSD,SP*b*mb*relEI/S0,SP*b*mb**2*relEI/S0,0,SP*b*mb*relEI/S0,SP*b*mb**2*relEI/S0,0,SP*b*mb*relEI/S0,SP*b*mb**2*relEI/S0,0,SP*b*mb*relEI/S0,SP*b*mb**2*relEI/S0,0,SP*b*mb*relA/S0,SP*b*mb**2*relA/S0,0,SP*b*mb/S0,SP*b*mb**2/S0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,SP*b*mb*relEI/S0,SP*b*mb*relA/S0,SP*b*mb/S0,0,0,0,],[0,0,kL,0,-kL - kQ - kfSD,krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,kL,kfSD,-kL - kQ - krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,kQ,kQ,-kL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,kL,0,0,-kL - kQ - kfSD,krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,kL,0,kfSD,-kL - kQ - krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,kL,kQ,kQ,-kL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,kL,0,0,-kL - kQ - kfSD,krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,kL,0,kfSD,-kL - kQ - krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,kL,kQ,kQ,-kL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,kL,0,0,-kL - kQ - kfSD,krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,kL,0,kfSD,-kL - kQ - krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,kL,kQ,kQ,-kL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,fA*kL,0,0,-gA - kQ - kfSD,krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,fA*kL,0,kfSD,-gA - kQ - krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,fA*kL,kQ,kQ,-gA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,kL*(1 - fA),0,0,0,0,0,-jQ - kH - kQ - kR - kfSD,krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,kL*(1 - fA),0,0,0,0,kfSD,-jQ - kH - kQ - kR - krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kL*(1 - fA),0,0,0,jQ + kQ,jQ + kQ,-kH - kR,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kH,kH,kH,-kHD - kHR,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + gA,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,kR,0,0,kHR,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP) - kfSD,RM*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + gA,gA,0,kR,kR,0,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + kfSD,RP*S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP) - krSD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kHD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kHD,0,0,],[-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP),-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP),-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*relEI/S0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*mb*relEI/S0,0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*relEI/S0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*mb*relEI/S0,0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*relEI/S0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*mb*relEI/S0,0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*relEI/S0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*mb*relEI/S0,0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*relA/S0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 - V1*b*mb*relA/S0,0,-V1*b/S0,-V1*b*mb/S0,0,0,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(SM + SP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,-kV - b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,0,0,0,0,0,0,0,0,-V1*b*relEI/S0,-V1*b*relA/S0,-V1*b/S0,0,0,0,],[0,0,0,0,-V2*b*relEI/S0,-V2*b*mb*relEI/S0,0,-V2*b*relEI/S0,-V2*b*mb*relEI/S0,0,-V2*b*relEI/S0,-V2*b*mb*relEI/S0,0,-V2*b*relEI/S0,-V2*b*mb*relEI/S0,0,-V2*b*relA/S0,-V2*b*mb*relA/S0,0,-V2*b/S0,-V2*b*mb/S0,0,0,0,0,0,kV,-kV - b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,0,0,0,0,0,0,0,-V2*b*relEI/S0,-V2*b*relA/S0,-V2*b/S0,0,0,0,],[0,0,0,0,-V3*b*relEI/S0,-V3*b*mb*relEI/S0,0,-V3*b*relEI/S0,-V3*b*mb*relEI/S0,0,-V3*b*relEI/S0,-V3*b*mb*relEI/S0,0,-V3*b*relEI/S0,-V3*b*mb*relEI/S0,0,-V3*b*relA/S0,-V3*b*mb*relA/S0,0,-V3*b/S0,-V3*b*mb/S0,0,0,0,0,0,0,kV,-kV - b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,0,0,0,0,0,0,-V3*b*relEI/S0,-V3*b*relA/S0,-V3*b/S0,0,0,0,],[0,0,0,0,-V4*b*relEI/S0,-V4*b*mb*relEI/S0,0,-V4*b*relEI/S0,-V4*b*mb*relEI/S0,0,-V4*b*relEI/S0,-V4*b*mb*relEI/S0,0,-V4*b*relEI/S0,-V4*b*mb*relEI/S0,0,-V4*b*relA/S0,-V4*b*mb*relA/S0,0,-V4*b/S0,-V4*b*mb/S0,0,0,0,0,0,0,0,kV,-kV - b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,0,0,0,0,0,-V4*b*relEI/S0,-V4*b*relA/S0,-V4*b/S0,0,0,0,],[0,0,0,0,-V5*b*relEI/S0,-V5*b*mb*relEI/S0,0,-V5*b*relEI/S0,-V5*b*mb*relEI/S0,0,-V5*b*relEI/S0,-V5*b*mb*relEI/S0,0,-V5*b*relEI/S0,-V5*b*mb*relEI/S0,0,-V5*b*relA/S0,-V5*b*mb*relA/S0,0,-V5*b/S0,-V5*b*mb/S0,0,0,0,0,0,0,0,0,kV,-kV - b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,0,0,0,0,-V5*b*relEI/S0,-V5*b*relA/S0,-V5*b/S0,0,0,0,],[0,0,0,0,-V6*b*relEI/S0,-V6*b*mb*relEI/S0,0,-V6*b*relEI/S0,-V6*b*mb*relEI/S0,0,-V6*b*relEI/S0,-V6*b*mb*relEI/S0,0,-V6*b*relEI/S0,-V6*b*mb*relEI/S0,0,-V6*b*relA/S0,-V6*b*mb*relA/S0,0,-V6*b/S0,-V6*b*mb/S0,0,0,0,0,0,0,0,0,0,kV,-kV - b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,0,0,0,-V6*b*relEI/S0,-V6*b*relA/S0,-V6*b/S0,0,0,0,],[0,0,0,0,-SV1*b*relEI/S0,-SV1*b*mb*relEI/S0,0,-SV1*b*relEI/S0,-SV1*b*mb*relEI/S0,0,-SV1*b*relEI/S0,-SV1*b*mb*relEI/S0,0,-SV1*b*relEI/S0,-SV1*b*mb*relEI/S0,0,-SV1*b*relA/S0,-SV1*b*mb*relA/S0,0,-SV1*b/S0,-SV1*b*mb/S0,0,0,0,0,0,0,0,0,0,0,kV*(1 - f0),-b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,0,0,-SV1*b*relEI/S0,-SV1*b*relA/S0,-SV1*b/S0,0,0,0,],[0,0,0,0,-SV2*b*hmutation1*relEI/S0,-SV2*b*hmutation1*mb*relEI/S0,0,-SV2*b*hmutation1*relEI/S0,-SV2*b*hmutation1*mb*relEI/S0,0,-SV2*b*hmutation1*relEI/S0,-SV2*b*hmutation1*mb*relEI/S0,0,-SV2*b*hmutation1*relEI/S0,-SV2*b*hmutation1*mb*relEI/S0,0,-SV2*b*hmutation1*relA/S0,-SV2*b*hmutation1*mb*relA/S0,0,-SV2*b*hmutation1/S0,-SV2*b*hmutation1*mb/S0,0,0,0,0,0,0,0,0,0,0,kV*(f0 - f1),0,-b*hmutation1*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,0,-SV2*b*hmutation1*relEI/S0,-SV2*b*hmutation1*relA/S0,-SV2*b*hmutation1/S0,0,0,0,],[0,0,0,0,-SV3*b*hmutation2*relEI/S0,-SV3*b*hmutation2*mb*relEI/S0,0,-SV3*b*hmutation2*relEI/S0,-SV3*b*hmutation2*mb*relEI/S0,0,-SV3*b*hmutation2*relEI/S0,-SV3*b*hmutation2*mb*relEI/S0,0,-SV3*b*hmutation2*relEI/S0,-SV3*b*hmutation2*mb*relEI/S0,0,-SV3*b*hmutation2*relA/S0,-SV3*b*hmutation2*mb*relA/S0,0,-SV3*b*hmutation2/S0,-SV3*b*hmutation2*mb/S0,0,0,0,0,0,0,0,0,0,0,kV*(f1 - f2),0,0,-b*hmutation2*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,-SV3*b*hmutation2*relEI/S0,-SV3*b*hmutation2*relA/S0,-SV3*b*hmutation2/S0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,f2*kV,0,0,0,0,0,0,0,0,0,0,],[0,0,0,0,b*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,b*mb*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,0,b*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,b*mb*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,0,b*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,b*mb*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,0,b*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,b*mb*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,0,b*relA*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,b*mb*relA*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,0,b*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,b*mb*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,0,0,0,0,0,b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,b*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,b*hmutation1*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,b*hmutation2*(IM + IV + mb*(AP*relA + IP + relEI*(E2P + E3P + E4P + E5P)) + relA*(AM + AV) + relEI*(E2M + E3M + E4M + E5M + EV))/S0,0,-kL/5 + b*relEI*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,b*relA*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,b*(SV1 + SV2*hmutation1 + SV3*hmutation2 + V1 + V2 + V3 + V4 + V5 + V6)/S0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,fA*kL/5,-gA,0,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kL*(1 - fA)/5,0,-kH - kR,0,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kH/25,-kHD - kHR,0,0,],[-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2,0,0,0,0,0,-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP),-S0*mu*(RM + RP)/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP)**2 + S0*mu/(AM + AP + E1M + E1P + E2M + E2P + E3M + E3P + E4M + E4P + E5M + E5P + RM + RP + SM + SP),0,0,0,0,0,0,0,0,0,0,0,0,gA,kR25,kHR,0,0,],[0,0,0,0,0,0,0,0,0,0,0,0,0,kL*(1 - fA),kL*(1 - fA),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kL*(1 - fA)/5,0,0,0,0,0,],])

def process(tSpan,par):
    sol = solve_ivp(fun=lambda t,z: RHS(t,z,par), t_span = (tSpan[0],tSpan[-1]), y0 = state, t_eval = tSpan, method='RK45',rtol=1e-4,atol=1e-4)
    newCases = par[nStages*3+2*nMutations+2]*(sol.y[-1,1:]-sol.y[-1,:-1])
    return newCases

def penalty(par):
    incidents = NYTNewCases
    tSpan = linspace(0, len(incidents), len(incidents)+1)
    output = process(tSpan, par)
    r = par[nStages*3+2*nMutations+3]
    prob = r/(r+output)
    prob[prob<=1e-10] = 1e-10
    prob[prob>=1-1e-10] = 1-1e-10
    v1 = prob[incidents>=0]
    v2 = incidents[incidents>=0]
    var = loggamma(v2+r)-loggamma(v2+1)-loggamma(r)+r*log(v1)+v2*log(1-v1)
    logLval = sum(var)
    return logLval
    
temp1 = load('/Users/amallela/Documents/bmab/inferences/NYC0-n4-parLog.npy')
temp2 = load('/Users/amallela/Documents/bmab/inferences/NYC0-n4-logLLog.npy')
par = temp1[where(temp2==amax(temp2))[0][0]]
dateIndex = ['2020-01-21', '2020-01-22', '2020-01-23', '2020-01-24', '2020-01-25', '2020-01-26', '2020-01-27', '2020-01-28', '2020-01-29', '2020-01-30', '2020-01-31', '2020-02-01', '2020-02-02', '2020-02-03', '2020-02-04', '2020-02-05', '2020-02-06', '2020-02-07', '2020-02-08', '2020-02-09', '2020-02-10', '2020-02-11', '2020-02-12', '2020-02-13', '2020-02-14', '2020-02-15', '2020-02-16', '2020-02-17', '2020-02-18', '2020-02-19', '2020-02-20', '2020-02-21', '2020-02-22', '2020-02-23', '2020-02-24', '2020-02-25', '2020-02-26', '2020-02-27', '2020-02-28', '2020-02-29', '2020-03-01', '2020-03-02', '2020-03-03', '2020-03-04', '2020-03-05', '2020-03-06', '2020-03-07', '2020-03-08', '2020-03-09', '2020-03-10', '2020-03-11', '2020-03-12', '2020-03-13', '2020-03-14', '2020-03-15', '2020-03-16', '2020-03-17', '2020-03-18', '2020-03-19', '2020-03-20', '2020-03-21', '2020-03-22', '2020-03-23', '2020-03-24', '2020-03-25', '2020-03-26', '2020-03-27', '2020-03-28', '2020-03-29', '2020-03-30', '2020-03-31', '2020-04-01', '2020-04-02', '2020-04-03', '2020-04-04', '2020-04-05', '2020-04-06', '2020-04-07', '2020-04-08', '2020-04-09', '2020-04-10', '2020-04-11', '2020-04-12', '2020-04-13', '2020-04-14', '2020-04-15', '2020-04-16', '2020-04-17', '2020-04-18', '2020-04-19', '2020-04-20', '2020-04-21', '2020-04-22', '2020-04-23', '2020-04-24', '2020-04-25', '2020-04-26', '2020-04-27', '2020-04-28', '2020-04-29', '2020-04-30', '2020-05-01', '2020-05-02', '2020-05-03', '2020-05-04', '2020-05-05', '2020-05-06', '2020-05-07', '2020-05-08', '2020-05-09', '2020-05-10', '2020-05-11', '2020-05-12', '2020-05-13', '2020-05-14', '2020-05-15', '2020-05-16', '2020-05-17', '2020-05-18', '2020-05-19', '2020-05-20', '2020-05-21', '2020-05-22', '2020-05-23', '2020-05-24', '2020-05-25', '2020-05-26', '2020-05-27', '2020-05-28', '2020-05-29', '2020-05-30', '2020-05-31', '2020-06-01', '2020-06-02', '2020-06-03', '2020-06-04', '2020-06-05', '2020-06-06', '2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', '2020-06-12', '2020-06-13', '2020-06-14', '2020-06-15', '2020-06-16', '2020-06-17', '2020-06-18', '2020-06-19', '2020-06-20', '2020-06-21', '2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25', '2020-06-26', '2020-06-27', '2020-06-28', '2020-06-29', '2020-06-30', '2020-07-01', '2020-07-02', '2020-07-03', '2020-07-04', '2020-07-05', '2020-07-06', '2020-07-07', '2020-07-08', '2020-07-09', '2020-07-10', '2020-07-11', '2020-07-12', '2020-07-13', '2020-07-14', '2020-07-15', '2020-07-16', '2020-07-17', '2020-07-18', '2020-07-19', '2020-07-20', '2020-07-21', '2020-07-22', '2020-07-23', '2020-07-24', '2020-07-25', '2020-07-26', '2020-07-27', '2020-07-28', '2020-07-29', '2020-07-30', '2020-07-31', '2020-08-01', '2020-08-02', '2020-08-03', '2020-08-04', '2020-08-05', '2020-08-06', '2020-08-07', '2020-08-08', '2020-08-09', '2020-08-10', '2020-08-11', '2020-08-12', '2020-08-13', '2020-08-14', '2020-08-15', '2020-08-16', '2020-08-17', '2020-08-18', '2020-08-19', '2020-08-20', '2020-08-21', '2020-08-22', '2020-08-23', '2020-08-24', '2020-08-25', '2020-08-26', '2020-08-27', '2020-08-28', '2020-08-29', '2020-08-30', '2020-08-31', '2020-09-01', '2020-09-02', '2020-09-03', '2020-09-04', '2020-09-05', '2020-09-06', '2020-09-07', '2020-09-08', '2020-09-09', '2020-09-10', '2020-09-11', '2020-09-12', '2020-09-13', '2020-09-14', '2020-09-15', '2020-09-16', '2020-09-17', '2020-09-18', '2020-09-19', '2020-09-20', '2020-09-21', '2020-09-22', '2020-09-23', '2020-09-24', '2020-09-25', '2020-09-26', '2020-09-27', '2020-09-28', '2020-09-29', '2020-09-30', '2020-10-01', '2020-10-02', '2020-10-03', '2020-10-04', '2020-10-05', '2020-10-06', '2020-10-07', '2020-10-08', '2020-10-09', '2020-10-10', '2020-10-11', '2020-10-12', '2020-10-13', '2020-10-14', '2020-10-15', '2020-10-16', '2020-10-17', '2020-10-18', '2020-10-19', '2020-10-20', '2020-10-21', '2020-10-22', '2020-10-23', '2020-10-24', '2020-10-25', '2020-10-26', '2020-10-27', '2020-10-28', '2020-10-29', '2020-10-30', '2020-10-31', '2020-11-01', '2020-11-02', '2020-11-03', '2020-11-04', '2020-11-05', '2020-11-06', '2020-11-07', '2020-11-08', '2020-11-09', '2020-11-10', '2020-11-11', '2020-11-12', '2020-11-13', '2020-11-14', '2020-11-15', '2020-11-16', '2020-11-17', '2020-11-18', '2020-11-19', '2020-11-20', '2020-11-21', '2020-11-22', '2020-11-23', '2020-11-24', '2020-11-25', '2020-11-26', '2020-11-27', '2020-11-28', '2020-11-29', '2020-11-30','2020-12-01', '2020-12-02', '2020-12-03', '2020-12-04', '2020-12-05', '2020-12-06', '2020-12-07', '2020-12-08', '2020-12-09', '2020-12-10', '2020-12-11', '2020-12-12', '2020-12-13', '2020-12-14', '2020-12-15', '2020-12-16', '2020-12-17', '2020-12-18', '2020-12-19', '2020-12-20', '2020-12-21', '2020-12-22', '2020-12-23', '2020-12-24', '2020-12-25', '2020-12-26', '2020-12-27', '2020-12-28', '2020-12-29', '2020-12-30','2020-12-31', '2021-01-01', '2021-01-02', '2021-01-03', '2021-01-04', '2021-01-05', '2021-01-06', '2021-01-07', '2021-01-08', '2021-01-09', '2021-01-10', '2021-01-11', '2021-01-12', '2021-01-13', '2021-01-14', '2021-01-15', '2021-01-16', '2021-01-17', '2021-01-18', '2021-01-19', '2021-01-20', '2021-01-21', '2021-01-22', '2021-01-23', '2021-01-24', '2021-01-25', '2021-01-26', '2021-01-27', '2021-01-28', '2021-01-29', '2021-01-30', '2021-01-31', '2021-02-01', '2021-02-02', '2021-02-03', '2021-02-04', '2021-02-05', '2021-02-06', '2021-02-07', '2021-02-08', '2021-02-09', '2021-02-10', '2021-02-11', '2021-02-12', '2021-02-13', '2021-02-14', '2021-02-15', '2021-02-16', '2021-02-17', '2021-02-18', '2021-02-19', '2021-02-20', '2021-02-21', '2021-02-22', '2021-02-23', '2021-02-24', '2021-02-25', '2021-02-26', '2021-02-27', '2021-02-28', '2021-02-29', '2021-03-01', '2021-03-02', '2021-03-03', '2021-03-04', '2021-03-05', '2021-03-06', '2021-03-07', '2021-03-08', '2021-03-09', '2021-03-10', '2021-03-11', '2021-03-12', '2021-03-13', '2021-03-14', '2021-03-15', '2021-03-16', '2021-03-17', '2021-03-18', '2021-03-19', '2021-03-20', '2021-03-21', '2021-03-22', '2021-03-23', '2021-03-24', '2021-03-25', '2021-03-26', '2021-03-27', '2021-03-28', '2021-03-29', '2021-03-30', '2021-03-31', '2021-04-01', '2021-04-02', '2021-04-03', '2021-04-04', '2021-04-05', '2021-04-06', '2021-04-07', '2021-04-08', '2021-04-09', '2021-04-10', '2021-04-11', '2021-04-12', '2021-04-13', '2021-04-14', '2021-04-15', '2021-04-16', '2021-04-17', '2021-04-18', '2021-04-19', '2021-04-20', '2021-04-21', '2021-04-22', '2021-04-23', '2021-04-24', '2021-04-25', '2021-04-26', '2021-04-27', '2021-04-28', '2021-04-29', '2021-04-30', '2021-05-01', '2021-05-02', '2021-05-03', '2021-05-04', '2021-05-05', '2021-05-06', '2021-05-07', '2021-05-08', '2021-05-09', '2021-05-10', '2021-05-11', '2021-05-12', '2021-05-13', '2021-05-14', '2021-05-15', '2021-05-16', '2021-05-17', '2021-05-18', '2021-05-19', '2021-05-20', '2021-05-21', '2021-05-22', '2021-05-23', '2021-05-24', '2021-05-25', '2021-05-26', '2021-05-27', '2021-05-28', '2021-05-29', '2021-05-30', '2021-05-31', '2021-06-01', '2021-06-02', '2021-06-03', '2021-06-04', '2021-06-05', '2021-06-06', '2021-06-07', '2021-06-08', '2021-06-09', '2021-06-10', '2021-06-11', '2021-06-12', '2021-06-13', '2021-06-14', '2021-06-15', '2021-06-16', '2021-06-17', '2021-06-18', '2021-06-19', '2021-06-20', '2021-06-21', '2021-06-22', '2021-06-23', '2021-06-24', '2021-06-25', '2021-06-26', '2021-06-27', '2021-06-28', '2021-06-29', '2021-06-30', '2021-07-01', '2021-07-02', '2021-07-03', '2021-07-04', '2021-07-05', '2021-07-06', '2021-07-07', '2021-07-08', '2021-07-09', '2021-07-10', '2021-07-11', '2021-07-12', '2021-07-13', '2021-07-14', '2021-07-15', '2021-07-16', '2021-07-17', '2021-07-18', '2021-07-19', '2021-07-20', '2021-07-21', '2021-07-22', '2021-07-23', '2021-07-24', '2021-07-25', '2021-07-26', '2021-07-27', '2021-07-28', '2021-07-29', '2021-07-30', '2021-07-31', '2021-08-01', '2021-08-02', '2021-08-03', '2021-08-04', '2021-08-05', '2021-08-06', '2021-08-07', '2021-08-08', '2021-08-09', '2021-08-10', '2021-08-11', '2021-08-12', '2021-08-13', '2021-08-14', '2021-08-15', '2021-08-16', '2021-08-17', '2021-08-18', '2021-08-19', '2021-08-20', '2021-08-21', '2021-08-22', '2021-08-23', '2021-08-24', '2021-08-25', '2021-08-26', '2021-08-27', '2021-08-28', '2021-08-29', '2021-08-30', '2021-08-31', '2021-09-01', '2021-09-02', '2021-09-03', '2021-09-04', '2021-09-05', '2021-09-06', '2021-09-07', '2021-09-08', '2021-09-09', '2021-09-10', '2021-09-11', '2021-09-12', '2021-09-13', '2021-09-14', '2021-09-15', '2021-09-16', '2021-09-17', '2021-09-18', '2021-09-19', '2021-09-20', '2021-09-21', '2021-09-22', '2021-09-23', '2021-09-24', '2021-09-25', '2021-09-26', '2021-09-27', '2021-09-28', '2021-09-29', '2021-09-30', '2021-10-01', '2021-10-02', '2021-10-03', '2021-10-04', '2021-10-05', '2021-10-06', '2021-10-07', '2021-10-08', '2021-10-09', '2021-10-10', '2021-10-11', '2021-10-12', '2021-10-13', '2021-10-14', '2021-10-15', '2021-10-16', '2021-10-17', '2021-10-18', '2021-10-19', '2021-10-20', '2021-10-21', '2021-10-22', '2021-10-23', '2021-10-24', '2021-10-25', '2021-10-26', '2021-10-27', '2021-10-28', '2021-10-29', '2021-10-30', '2021-10-31']
fig = plt.figure(figsize=(24,24))
ax = plt.subplot(1,1,1)
ax.plot(tSpan,process(tSpanSimulation,par), linewidth=20.0, zorder = 1)    
ax.scatter(tSpan, NYTNewCases, 500, marker='+', linewidth = 4.0, color = 'black', zorder = 2)
ax.set_xticks(range(0,648,15))
_ = ax.set_xticklabels(dateIndex[::15], rotation=90)
ax.set_ylim(0,55000)
ax.set_title('NYC', fontsize = 30)
ax.set_ylabel('Confirmed case counts')
fig.tight_layout()
plt.savefig('NYC-mle-plot.pdf',dpi=1200)
plt.close()