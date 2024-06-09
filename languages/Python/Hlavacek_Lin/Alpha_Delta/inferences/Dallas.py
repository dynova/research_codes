from numpy import *
from scipy import stats
from scipy.integrate import odeint
from scipy.special import loggamma
from datetime import datetime
from sys import argv
from numba import njit

nStages = 5
nMutations = 2
locationIndex = 3
WarmStart = False
nms = ['NYC','LA','Chicago','Dallas','Houston','DC','Philadelphia','Miami','Atlanta','Boston','Phoenix','SanFrancisco','Riverside','Detroit','Seattle']
population = [19216182,13214799,9458539,7573136,7066141,6280487,6166488,6102434,6020364,4948203,4873019,4731803,4650631,4319629,3979845]
prefix = nms[locationIndex]
totalPopulation = population[locationIndex]
filePrefix = prefix+str(int(argv[1]))+'-n'+str(nStages)
NYTNewCases = array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,5,1,6,3,6,5,1,12,20,46,43,53,52,50,21,64,191,183,97,55,163,135,277,196,198,160,134,170,181,275,228,207,196,206,199,174,193,146,380,203,129,125,243,250,260,264,283,238,204,249,251,316,404,370,325,317,391,419,498,461,377,761,368,374,565,451,323,336,332,390,407,341,353,344,257,231,184,348,409,401,371,373,315,281,388,443,479,551,440,361,565,543,598,626,759,562,420,641,607,834,859,799,778,948,738,958,1487,1045,1006,1365,1118,1088,1677,1551,1515,1597,1482,1818,2342,2363,2325,2229,2239,1949,1787,1735,2307,2502,1957,2268,1648,1657,1634,2082,1823,1616,1442,2236,1269,913,2259,1961,1589,3064,67,2226,863,2437,2867,1424,1492,1840,1737,1274,1582,1332,1579,2871,2866,6072,2399,2379,920,878,1366,1543,640,491,1018,1183,771,812,1271,281,1201,1339,815,439,944,971,416,115,950,1204,494,909,675,329,1389,1189,998,1081,1113,718,1319,812,989,1072,1655,482,721,1996,857,1729,874,629,492,2983,1275,847,883,2449,1099,1062,1163,1294,879,1529,1550,1625,1516,1817,1168,936,1265,1822,1315,1920,2337,1524,1296,1661,1741,1906,1097,2625,1752,2440,2095,2168,1675,3724,1258,2535,3321,3442,2994,2353,3214,3396,3256,3069,2063,4016,5000,4296,4421,4124,2552,6226,3773,897,626,1609,5506,4834,3742,6499,6159,3539,4274,4079,4498,4324,4432,5016,12974,4127,1900,5526,6913,6342,6792,8117,5539,1385,5927,7355,7241,1930,74,2920,7506,10318,6419,6030,7251,1247,4071,8154,9208,10210,8892,7553,7547,5919,6108,8359,9861,8180,7303,6874,4802,5068,5235,5461,10277,5689,6286,3873,3574,6152,4750,6624,6475,5240,2774,3814,6807,5845,4982,4475,4453,2547,2590,3869,4155,4277,3889,3176,1904,1806,890,1735,2429,878,1170,1178,351,1916,2571,2470,2556,2594,2078,189,2557,1925,2251,2766,1564,948,425,1207,2010,1888,1510,1367,488,342,1513,1236,907,1077,883,356,142,1224,994,694,652,1041,506,172,1287,818,1108,919,444,331,139,654,1421,796,636,915,466,278,671,1383,640,752,718,456,353,601,1247,654,740,641,322,631,430,1251,953,879,861,599,506,665,883,569,702,568,337,454,439,715,287,702,663,204,268,402,1216,340,408,290,237,291,336,596,372,243,501,123,51,50,51,521,428,336,110,175,229,208,450,356,254,941,156,254,474,371,525,301,271,286,200,402,283,480,389,91,63,334,371,451,193,559,331,37,27,387,557,600,641,363,142,446,763,1001,1052,1133,519,966,934,1145,1396,2006,1291,767,777,1274,2107,2788,2610,2227,1377,1186,1934,2789,3976,2866,3396,1405,886,2700,3775,5634,3443,3124,1923,1967,3101,4920,7079,3415,3060,2602,1257,2463,4288,6580,3942,5097,1787,749,4827,4526,7408,4882,6348,2972,236,898,3745,9247,10131,4726,2500,459,8590,6033,7207,4531,4560,2948,191,4634,3384,3462,2247,5070,1234,3683,3748,3465,4741,2774,3088,1639,202,3030,3165,3964,2802,2254,1422,80,2579,2866,3026,2082,2201,1337,104,1907,2012,2284,1876,1456,1044,141,1559,1881,1486,2193,1310,872,105])
vr_daily=load('vac.npy')
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
ts = nan*zeros(nStages+1)                                                

ts[0]= 1.54644212e+01            
ts[1]= 5.21444582e+01
ts[2]= 8.40664694e+01
ts[3]= 2.58820063e+01
ts[4]= 9.99140809e+01
ts[5]= 6.21781596e+00

tmutation = nan*zeros(nMutations)

mutatinMultiplier = nan*zeros(nMutations)

lambs = nan*zeros(nStages)
lambs[0] = 9.47018170e+00         
lambs[1] = 1.76433900e+00
lambs[2] = 2.03510847e-01
lambs[3] = 4.57584442e+00
lambs[4] = 1.45256022e+00

fPs = nan*zeros(nStages)
fPs[0] = 4.96637848e-01         
fPs[1] = 4.20208777e-01
fPs[2] = 4.21891837e-01
fPs[3] = 1.98647946e-01
fPs[4] = 2.30659542e-01

tmutation[0] = 3.83530833e+02   
tmutation[1] = 1.28087762e+02   
mutatinMultiplier[0] = 1.42439918e+00   
mutatinMultiplier[1] = 2.11290739e+00   

b = 3.87259496e-01   
fDI = 9.97760129e-01   
r = 9.09420557e+00

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
parRest = [b, fDI, r, S0, mb, relEI, relA, kL, gI, fA, fH, kH, kR, kQ, jQ, gA, gH, kHD, kHR, kV, kR25, f0, f1, f2]
par = []

for i in range(nStages+1):
    par.append(ts[i])
    
for i in range(nStages):
    par.append(lambs[i])
    
for i in range(nStages):
    par.append(fPs[i])

for i in range(nMutations):
    par.append(tmutation[i])
    
for i in range(nMutations):
    par.append(mutatinMultiplier[i])    

for i in range(len(parRest)):
    par.append(parRest[i])
    
par = array(par)
    
freeParameterInd = range(4+nStages*3+nMutations*2)

parameterBounds = []

for i in range(nStages+1):
    parameterBounds.append([0.0, inf])
    
for i in range(nStages):
    parameterBounds.append([0.0, 10.0])
    
for i in range(nStages):
    parameterBounds.append([0.0, 1.0])

for i in range(nMutations):
    parameterBounds.append([0.0, inf])
    
for i in range(nMutations):
    parameterBounds.append([1.0, 20])
    
parameterBounds.append([0.0, inf])
parameterBounds.append([0.0, 1.0])
parameterBounds.append([0.0, inf])

parameterBounds = array(parameterBounds)

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

@njit
def RHS(state,t,par):

    SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,V1,V2,V3,V4,V5,V6,SV1,SV2,SV3,SV4,EV,AV,IV,HV,RV,cI= state
    
    ts = par[:nStages+1]
    lambs = par[nStages+1:2*nStages+1]
    fPs = par[2*nStages+1:3*nStages+1]
    tmutation = par[3*nStages+1:3*nStages+1+nMutations]
    mutationMultiplier = par[3*nStages+1+nMutations:3*nStages+1+2*nMutations]
    
    b, fDI, r, S0, mb, relEI, relA, kL, gI, fA, fH, kH, kR, kQ, jQ, gA, gH, kHD, kHR, kV, kR25, f0, f1, f2 = par[nStages*3+2*nMutations+1:]
    
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
def Jacobian(state,t,par):

    SM,SP,E1M,E1P,E2M,E2P,E2Q,E3M,E3P,E3Q,E4M,E4P,E4Q,E5M,E5P,E5Q,AM,AP,AQ,IM,IP,IQ,IH,RM,RP,D,V1,V2,V3,V4,V5,V6,SV1,SV2,SV3,SV4,EV,AV,IV,HV,RV,cI= state
    
    ts = par[:nStages+1]
    lambs = par[nStages+1:2*nStages+1]
    fPs = par[2*nStages+1:3*nStages+1]
    tmutation = par[3*nStages+1:3*nStages+1+nMutations]
    mutationMultiplier = par[3*nStages+1+nMutations:3*nStages+1+2*nMutations]
    
    b, fDI, r, S0, mb, relEI, relA, kL, gI, fA, fH, kH, kR, kQ, jQ, gA, gH, kHD, kHR, kV, kR25, f0, f1, f2 = par[nStages*3+2*nMutations+1:]
    
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
    sol = odeint(RHS,state,tSpan,args=(par,),Dfun=Jacobian,rtol=1e-4,atol=1e-4)
    newCases = par[nStages*3+2*nMutations+2]*(sol[1:,-1]-sol[:-1,-1])
    return newCases

def logL(par, incidents):
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

def generateBinomialNoise(timeseries, par):
    output = copy(timeseries)
    r = par[nStages*3+2*nMutations+3]
    prob = r/(r+timeseries)
    prob[prob<=1e-10] = 1e-10
    prob[prob>=1-1e-10] = 1-1e-10
    output = stats.nbinom.rvs(n=r, p=prob, size=len(timeseries))
    return output

def MCMC(iniPar, freeParameterInd, parameterBounds, diffMatrix, diff, maxIter, iterSaving, incidents, prefix, cutoffIter=150000, iterBurning=50000, iterTurningOnAdaptive=50000):

    totalParN = len(freeParameterInd)
    stabilizingCov = 0.001*eye(totalParN)

    if len(freeParameterInd)!= len(parameterBounds):
        print('The dimension of the bounds is not consistent with the free parameter list.')
        return

    tSpan = linspace(0, len(incidents)+14, len(incidents)+15)  # forecast for additional 2 weeks

    allTrajectoryLog = zeros((maxIter, len(tSpan)-1))  # An array that stores all the trajectories

    qtlMark = linspace(0, 100, 21)
    qtlMark[0] = 2.5
    qtlMark[-1] = 97.5                          # The quantiles/percentiles we will be plotting

    qtlLog = zeros((len(qtlMark),len(tSpan)-1))   # An array that stores all the quantiles a differet times

    parLog = zeros((maxIter,totalParN ))        # An array thar stores the parameters

    acceptLog = zeros((maxIter,))               # An array that stores the acceptance

    logLLog = -1E30*ones((maxIter,))            # An array that stores the likelihood

    par = copy(iniPar)

    parLog[0,:] = par[freeParameterInd]

    oldLogL = logL(par, incidents)

    logLLog[0] = oldLogL

    maxLogL = copy(oldLogL)

    output = process(tSpan, par)

    allTrajectoryLog[0,:] = generateBinomialNoise(output, par)

    MLE = output[:]

    for index in range(1,maxIter):

        ifCut = True

        while ifCut:

            ifCut = False

            proposedPar = copy(par)

            proposedPar[freeParameterInd] = proposedPar[freeParameterInd] + diff*random.multivariate_normal(mean=zeros((totalParN,)), cov=diffMatrix)

            # start prior cutoff here

            for kk in range(totalParN):

                ifCut = ifCut | (proposedPar[freeParameterInd[kk]] < parameterBounds[kk][0])
                ifCut = ifCut | (proposedPar[freeParameterInd[kk]] > parameterBounds[kk][1])


        newLogL= logL(proposedPar, incidents)

        accept = False

        if newLogL > oldLogL:

            accept = True
            alpha = 1

        else:

            alpha = exp((newLogL-oldLogL))

            if random.random() < alpha :

                accept = True

        if accept:

            tempOutput  = process(tSpan, proposedPar)

            if amin(tempOutput>=0):

                par = copy(proposedPar)
                oldLogL = copy(newLogL)
                acceptLog[index] = 1
                output = copy(tempOutput)

        logLLog[index] = oldLogL

        allTrajectoryLog[index,:] = generateBinomialNoise(output, par)

        if oldLogL>maxLogL:

            maxLogL = copy(oldLogL)
            MLE = output[:]

        parLog[index,:] = par[freeParameterInd]

        if index >= iterBurning+iterTurningOnAdaptive:

            if index == iterBurning+iterTurningOnAdaptive:

                mu = reshape(mean(parLog[iterBurning:index,:],axis=0), [1, totalParN])  # compute the mean parameters along the past chain
                diffMatrix = matmul(parLog[iterBurning:index,:].T, parLog[iterBurning:index,:])/(index-iterBurning)-matmul(mu.T, mu)+stabilizingCov
                diff = 2.38**2/totalParN

            mu = mu + (1./(1+index))*(parLog[index,:] - mu)
            diffVector = reshape(parLog[index,:]-mu, [1, totalParN])
            diffMatrix = diffMatrix + (1./(1+index))*(matmul(diffVector.T, diffVector)+stabilizingCov-diffMatrix)
            diff = exp( log(diff) + (1./(1+index-iterTurningOnAdaptive-iterBurning))*(alpha-0.234))
                
        if (index+1)%iterSaving==0:

            savetxt(prefix+'-epoch.txt',atleast_1d(index//iterSaving+1),fmt='%i')
            save(prefix+'-MLE.npy', MLE)
            save(prefix+'-logLLog.npy', logLLog)
            save(prefix+'-parLog.npy', parLog)
            save(prefix+'-acceptLog.npy', acceptLog)
            save(prefix+'-percentiles.npy', qtlLog)
            save(prefix+'-diffMatrix.npy', diffMatrix)
            save(prefix+'-diff.npy', diff)

            if (index+1)>cutoffIter:
            
                for i in range(len(tSpan)-1):

                    qtlLog[:,i] = percentile(allTrajectoryLog[cutoffIter:index+1,i], qtlMark)            
                            
    return logLLog, parLog, acceptLog, allTrajectoryLog        
        
if WarmStart:
    logLLog = load('Dallas-n5-logLLog.npy')
    parLog = load('Dallas-n5-parLog.npy')
    MLE = load('Dallas-n5-MLE.npy')
    diffMatrix = load('Dallas-n5-diffMatrix.npy')
    diff = load('Dallas-n5-diff.npy')
    index = where(logLLog==amax(logLLog))[0][0]
    parResume = copy(par)
    parResume[freeparameterInd] = parLog[index]
    par = copy(parResume)
    nbinomIteration = 1000000
    iterBurning = 1000
    iterTurningOnAdaptive = 1000
    cutoffIter = 2000
else:
    diffMatrix = diag(array(par)[freeParameterInd])
    diff = 2.38**2/len(freeParameterInd)
    nbinomIteration = 400000
    iterBurning = 50000
    iterTurningOnAdaptive = 50000
    cutoffIter = 150000

logLLog, parLog, acceptLog, allTrajectory = MCMC(par, freeParameterInd, parameterBounds, diffMatrix, diff, nbinomIteration, 1000, NYTNewCases, prefix=filePrefix, cutoffIter=cutoffIter, iterBurning=iterBurning, iterTurningOnAdaptive=iterTurningOnAdaptive)