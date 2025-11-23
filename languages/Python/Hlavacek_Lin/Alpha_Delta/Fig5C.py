from numpy import *
from scipy.integrate import odeint, solve_ivp
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
nStages = 5
nMutations = 2
locationIndex = 10
nms = ['NYC','LA','Chicago','Dallas','Houston','DC','Philadelphia','Miami','Atlanta','Boston','Phoenix','SanFrancisco','Riverside','Detroit','Seattle']
population = [19216182,13214799,9458539,7573136,7066141,6280487,6166488,6102434,6020364,4948203,4873019,4731803,4650631,4319629,3979845]
prefix = nms[locationIndex]
totalPopulation = population[locationIndex]
NYTNewCases = array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,5,1,6,3,6,5,1,12,20,46,43,53,52,50,21,64,191,183,97,55,163,135,277,196,198,160,134,170,181,275,228,207,196,206,199,174,193,146,380,203,129,125,243,250,260,264,283,238,204,249,251,316,404,370,325,317,391,419,498,461,377,761,368,374,565,451,323,336,332,390,407,341,353,344,257,231,184,348,409,401,371,373,315,281,388,443,479,551,440,361,565,543,598,626,759,562,420,641,607,834,859,799,778,948,738,958,1487,1045,1006,1365,1118,1088,1677,1551,1515,1597,1482,1818,2342,2363,2325,2229,2239,1949,1787,1735,2307,2502,1957,2268,1648,1657,1634,2082,1823,1616,1442,2236,1269,913,2259,1961,1589,3064,67,2226,863,2437,2867,1424,1492,1840,1737,1274,1582,1332,1579,2871,2866,6072,2399,2379,920,878,1366,1543,640,491,1018,1183,771,812,1271,281,1201,1339,815,439,944,971,416,115,950,1204,494,909,675,329,1389,1189,998,1081,1113,718,1319,812,989,1072,1655,482,721,1996,857,1729,874,629,492,2983,1275,847,883,2449,1099,1062,1163,1294,879,1529,1550,1625,1516,1817,1168,936,1265,1822,1315,1920,2337,1524,1296,1661,1741,1906,1097,2625,1752,2440,2095,2168,1675,3724,1258,2535,3321,3442,2994,2353,3214,3396,3256,3069,2063,4016,5000,4296,4421,4124,2552,6226,3773,897,626,1609,5506,4834,3742,6499,6159,3539,4274,4079,4498,4324,4432,5016,12974,4127,1900,5526,6913,6342,6792,8117,5539,1385,5927,7355,7241,1930,74,2920,7506,10318,6419,6030,7251,1247,4071,8154,9208,10210,8892,7553,7547,5919,6108,8359,9861,8180,7303,6874,4802,5068,5235,5461,10277,5689,6286,3873,3574,6152,4750,6624,6475,5240,2774,3814,6807,5845,4982,4475,4453,2547,2590,3869,4155,4277,3889,3176,1904,1806,890,1735,2429,878,1170,1178,351,1916,2571,2470,2556,2594,2078,189,2557,1925,2251,2766,1564,948,425,1207,2010,1888,1510,1367,488,342,1513,1236,907,1077,883,356,142,1224,994,694,652,1041,506,172,1287,818,1108,919,444,331,139,654,1421,796,636,915,466,278,671,1383,640,752,718,456,353,601,1247,654,740,641,322,631,430,1251,953,879,861,599,506,665,883,569,702,568,337,454,439,715,287,702,663,204,268,402,1216,340,408,290,237,291,336,596,372,243,501,123,51,50,51,521,428,336,110,175,229,208,450,356,254,941,156,254,474,371,525,301,271,286,200,402,283,480,389,91,63,334,371,451,193,559,331,37,27,387,557,600,641,363,142,446,763,1001,1052,1133,519,966,934,1145,1396,2006,1291,767,777,1274,2107,2788,2610,2227,1377,1186,1934,2789,3976,2866,3396,1405,886,2700,3775,5634,3443,3124,1923,1967,3101,4920,7079,3415,3060,2602,1257,2463,4288,6580,3942,5097,1787,749,4827,4526,7408,4882,6348,2972,236,898,3745,9247,10131,4726,2500,459,8590,6033,7207,4531,4560,2948,191,4634,3384,3462,2247,5070,1234,3683,3748,3465,4741,2774,3088,1639,202,3030,3165,3964,2802,2254,1422,80,2579,2866,3026,2082,2201,1337,104,1907,2012,2284,1876,1456,1044,141,1559,1881,1486,2193,1310,872,105])
vr_daily=load('vac.npy')
x1 = linspace(1, len(vr_daily[locationIndex])+2, len(vr_daily[locationIndex])+2) 
y1 = vr_daily[locationIndex]
y1 = append(y1, [0,0])

l = load('inferences/Phoenix29-n5-logLLog.npy')
x = load('inferences/Phoenix29-n5-parLog.npy')
bestFitPar = x[where(l==amax(l))[0][0]]

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

ts = nan*zeros(nStages+1)                                                                
lambs = nan*zeros(nStages)
fPs = nan*zeros(nStages)
tmutation = nan*zeros(nMutations)
mutatinMultiplier = nan*zeros(nMutations)
ts[0]= bestFitPar[0]         
ts[1]= bestFitPar[1]
ts[2]= bestFitPar[2]
ts[3]= bestFitPar[3]
ts[4]= bestFitPar[4]
ts[5]= bestFitPar[5]
lambs[0] = bestFitPar[6]        
lambs[1] = bestFitPar[7]
lambs[2] = bestFitPar[8]
lambs[3] = bestFitPar[9]
lambs[4] = bestFitPar[10]
fPs[0] = bestFitPar[11]
fPs[1] = bestFitPar[12]
fPs[2] = bestFitPar[13]
fPs[3] = bestFitPar[14]
fPs[4] = bestFitPar[15]
tmutation[0] = bestFitPar[16]   
tmutation[1] = bestFitPar[17]   
mutatinMultiplier[0] = bestFitPar[18]   
mutatinMultiplier[1] = bestFitPar[19]   
b = bestFitPar[20]   
fDI = bestFitPar[21]   
r = bestFitPar[22]

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
    return sol
    
mult = 1
tSpan = linspace(0, len(NYTNewCases)-1, len(NYTNewCases)*mult)
tmutation = bestFitPar[3*nStages+1:3*nStages+1+nMutations]
cumsumtmutation = cumsum(tmutation)[:]
cumsumtmutation = append(cumsumtmutation,inf)
hmutation1 = (cumsumtmutation[0]<tSpan)*1.0
hmutation2 = (cumsumtmutation[1]<tSpan)*1.0
sol = process(tSpan,par)     
susceptible_and_unvaccinated = sum(sol[:,:2],1)
susceptible_and_vaccinated = sum(sol[:,26:33],1) + hmutation1 * sol[:,33] + hmutation2 * sol[:,34]
actively_infected = sum(sol[:,2:23],1) + sum(sol[:,36:40],1)
removed_and_unvaccinated = sum(sol[:,23:26],1)
removed_and_vaccinated = (1 - hmutation1) * sol[:,33] + (1 - hmutation2) * sol[:,34] + sol[:,35] + sol[:,40]

plt.rcParams.update({'font.family':'sans','font.size':12})
dateIndex = ['2020-01-21', '2020-01-22', '2020-01-23', '2020-01-24', '2020-01-25', '2020-01-26', '2020-01-27', '2020-01-28', '2020-01-29', '2020-01-30', '2020-01-31', '2020-02-01', '2020-02-02', '2020-02-03', '2020-02-04', '2020-02-05', '2020-02-06', '2020-02-07', '2020-02-08', '2020-02-09', '2020-02-10', '2020-02-11', '2020-02-12', '2020-02-13', '2020-02-14', '2020-02-15', '2020-02-16', '2020-02-17', '2020-02-18', '2020-02-19', '2020-02-20', '2020-02-21', '2020-02-22', '2020-02-23', '2020-02-24', '2020-02-25', '2020-02-26', '2020-02-27', '2020-02-28', '2020-02-29', '2020-03-01', '2020-03-02', '2020-03-03', '2020-03-04', '2020-03-05', '2020-03-06', '2020-03-07', '2020-03-08', '2020-03-09', '2020-03-10', '2020-03-11', '2020-03-12', '2020-03-13', '2020-03-14', '2020-03-15', '2020-03-16', '2020-03-17', '2020-03-18', '2020-03-19', '2020-03-20', '2020-03-21', '2020-03-22', '2020-03-23', '2020-03-24', '2020-03-25', '2020-03-26', '2020-03-27', '2020-03-28', '2020-03-29', '2020-03-30', '2020-03-31', '2020-04-01', '2020-04-02', '2020-04-03', '2020-04-04', '2020-04-05', '2020-04-06', '2020-04-07', '2020-04-08', '2020-04-09', '2020-04-10', '2020-04-11', '2020-04-12', '2020-04-13', '2020-04-14', '2020-04-15', '2020-04-16', '2020-04-17', '2020-04-18', '2020-04-19', '2020-04-20', '2020-04-21', '2020-04-22', '2020-04-23', '2020-04-24', '2020-04-25', '2020-04-26', '2020-04-27', '2020-04-28', '2020-04-29', '2020-04-30', '2020-05-01', '2020-05-02', '2020-05-03', '2020-05-04', '2020-05-05', '2020-05-06', '2020-05-07', '2020-05-08', '2020-05-09', '2020-05-10', '2020-05-11', '2020-05-12', '2020-05-13', '2020-05-14', '2020-05-15', '2020-05-16', '2020-05-17', '2020-05-18', '2020-05-19', '2020-05-20', '2020-05-21', '2020-05-22', '2020-05-23', '2020-05-24', '2020-05-25', '2020-05-26', '2020-05-27', '2020-05-28', '2020-05-29', '2020-05-30', '2020-05-31', '2020-06-01', '2020-06-02', '2020-06-03', '2020-06-04', '2020-06-05', '2020-06-06', '2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', '2020-06-12', '2020-06-13', '2020-06-14', '2020-06-15', '2020-06-16', '2020-06-17', '2020-06-18', '2020-06-19', '2020-06-20', '2020-06-21', '2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25', '2020-06-26', '2020-06-27', '2020-06-28', '2020-06-29', '2020-06-30', '2020-07-01', '2020-07-02', '2020-07-03', '2020-07-04', '2020-07-05', '2020-07-06', '2020-07-07', '2020-07-08', '2020-07-09', '2020-07-10', '2020-07-11', '2020-07-12', '2020-07-13', '2020-07-14', '2020-07-15', '2020-07-16', '2020-07-17', '2020-07-18', '2020-07-19', '2020-07-20', '2020-07-21', '2020-07-22', '2020-07-23', '2020-07-24', '2020-07-25', '2020-07-26', '2020-07-27', '2020-07-28', '2020-07-29', '2020-07-30', '2020-07-31', '2020-08-01', '2020-08-02', '2020-08-03', '2020-08-04', '2020-08-05', '2020-08-06', '2020-08-07', '2020-08-08', '2020-08-09', '2020-08-10', '2020-08-11', '2020-08-12', '2020-08-13', '2020-08-14', '2020-08-15', '2020-08-16', '2020-08-17', '2020-08-18', '2020-08-19', '2020-08-20', '2020-08-21', '2020-08-22', '2020-08-23', '2020-08-24', '2020-08-25', '2020-08-26', '2020-08-27', '2020-08-28', '2020-08-29', '2020-08-30', '2020-08-31', '2020-09-01', '2020-09-02', '2020-09-03', '2020-09-04', '2020-09-05', '2020-09-06', '2020-09-07', '2020-09-08', '2020-09-09', '2020-09-10', '2020-09-11', '2020-09-12', '2020-09-13', '2020-09-14', '2020-09-15', '2020-09-16', '2020-09-17', '2020-09-18', '2020-09-19', '2020-09-20', '2020-09-21', '2020-09-22', '2020-09-23', '2020-09-24', '2020-09-25', '2020-09-26', '2020-09-27', '2020-09-28', '2020-09-29', '2020-09-30', '2020-10-01', '2020-10-02', '2020-10-03', '2020-10-04', '2020-10-05', '2020-10-06', '2020-10-07', '2020-10-08', '2020-10-09', '2020-10-10', '2020-10-11', '2020-10-12', '2020-10-13', '2020-10-14', '2020-10-15', '2020-10-16', '2020-10-17', '2020-10-18', '2020-10-19', '2020-10-20', '2020-10-21', '2020-10-22', '2020-10-23', '2020-10-24', '2020-10-25', '2020-10-26', '2020-10-27', '2020-10-28', '2020-10-29', '2020-10-30', '2020-10-31', '2020-11-01', '2020-11-02', '2020-11-03', '2020-11-04', '2020-11-05', '2020-11-06', '2020-11-07', '2020-11-08', '2020-11-09', '2020-11-10', '2020-11-11', '2020-11-12', '2020-11-13', '2020-11-14', '2020-11-15', '2020-11-16', '2020-11-17', '2020-11-18', '2020-11-19', '2020-11-20', '2020-11-21', '2020-11-22', '2020-11-23', '2020-11-24', '2020-11-25', '2020-11-26', '2020-11-27', '2020-11-28', '2020-11-29', '2020-11-30','2020-12-01', '2020-12-02', '2020-12-03', '2020-12-04', '2020-12-05', '2020-12-06', '2020-12-07', '2020-12-08', '2020-12-09', '2020-12-10', '2020-12-11', '2020-12-12', '2020-12-13', '2020-12-14', '2020-12-15', '2020-12-16', '2020-12-17', '2020-12-18', '2020-12-19', '2020-12-20', '2020-12-21', '2020-12-22', '2020-12-23', '2020-12-24', '2020-12-25', '2020-12-26', '2020-12-27', '2020-12-28', '2020-12-29', '2020-12-30','2020-12-31', '2021-01-01', '2021-01-02', '2021-01-03', '2021-01-04', '2021-01-05', '2021-01-06', '2021-01-07', '2021-01-08', '2021-01-09', '2021-01-10', '2021-01-11', '2021-01-12', '2021-01-13', '2021-01-14', '2021-01-15', '2021-01-16', '2021-01-17', '2021-01-18', '2021-01-19', '2021-01-20', '2021-01-21', '2021-01-22', '2021-01-23', '2021-01-24', '2021-01-25', '2021-01-26', '2021-01-27', '2021-01-28', '2021-01-29', '2021-01-30', '2021-01-31', '2021-02-01', '2021-02-02', '2021-02-03', '2021-02-04', '2021-02-05', '2021-02-06', '2021-02-07', '2021-02-08', '2021-02-09', '2021-02-10', '2021-02-11', '2021-02-12', '2021-02-13', '2021-02-14', '2021-02-15', '2021-02-16', '2021-02-17', '2021-02-18', '2021-02-19', '2021-02-20', '2021-02-21', '2021-02-22', '2021-02-23', '2021-02-24', '2021-02-25', '2021-02-26', '2021-02-27', '2021-02-28', '2021-02-29', '2021-03-01', '2021-03-02', '2021-03-03', '2021-03-04', '2021-03-05', '2021-03-06', '2021-03-07', '2021-03-08', '2021-03-09', '2021-03-10', '2021-03-11', '2021-03-12', '2021-03-13', '2021-03-14', '2021-03-15', '2021-03-16', '2021-03-17', '2021-03-18', '2021-03-19', '2021-03-20', '2021-03-21', '2021-03-22', '2021-03-23', '2021-03-24', '2021-03-25', '2021-03-26', '2021-03-27', '2021-03-28', '2021-03-29', '2021-03-30', '2021-03-31', '2021-04-01', '2021-04-02', '2021-04-03', '2021-04-04', '2021-04-05', '2021-04-06', '2021-04-07', '2021-04-08', '2021-04-09', '2021-04-10', '2021-04-11', '2021-04-12', '2021-04-13', '2021-04-14', '2021-04-15', '2021-04-16', '2021-04-17', '2021-04-18', '2021-04-19', '2021-04-20', '2021-04-21', '2021-04-22', '2021-04-23', '2021-04-24', '2021-04-25', '2021-04-26', '2021-04-27', '2021-04-28', '2021-04-29', '2021-04-30', '2021-05-01', '2021-05-02', '2021-05-03', '2021-05-04', '2021-05-05', '2021-05-06', '2021-05-07', '2021-05-08', '2021-05-09', '2021-05-10', '2021-05-11', '2021-05-12', '2021-05-13', '2021-05-14', '2021-05-15', '2021-05-16', '2021-05-17', '2021-05-18', '2021-05-19', '2021-05-20', '2021-05-21', '2021-05-22', '2021-05-23', '2021-05-24', '2021-05-25', '2021-05-26', '2021-05-27', '2021-05-28', '2021-05-29', '2021-05-30', '2021-05-31', '2021-06-01', '2021-06-02', '2021-06-03', '2021-06-04', '2021-06-05', '2021-06-06', '2021-06-07', '2021-06-08', '2021-06-09', '2021-06-10', '2021-06-11', '2021-06-12', '2021-06-13', '2021-06-14', '2021-06-15', '2021-06-16', '2021-06-17', '2021-06-18', '2021-06-19', '2021-06-20', '2021-06-21', '2021-06-22', '2021-06-23', '2021-06-24', '2021-06-25', '2021-06-26', '2021-06-27', '2021-06-28', '2021-06-29', '2021-06-30', '2021-07-01', '2021-07-02', '2021-07-03', '2021-07-04', '2021-07-05', '2021-07-06', '2021-07-07', '2021-07-08', '2021-07-09', '2021-07-10', '2021-07-11', '2021-07-12', '2021-07-13', '2021-07-14', '2021-07-15', '2021-07-16', '2021-07-17', '2021-07-18', '2021-07-19', '2021-07-20', '2021-07-21', '2021-07-22', '2021-07-23', '2021-07-24', '2021-07-25', '2021-07-26', '2021-07-27', '2021-07-28', '2021-07-29', '2021-07-30', '2021-07-31', '2021-08-01', '2021-08-02', '2021-08-03', '2021-08-04', '2021-08-05', '2021-08-06', '2021-08-07', '2021-08-08', '2021-08-09', '2021-08-10', '2021-08-11', '2021-08-12', '2021-08-13', '2021-08-14', '2021-08-15', '2021-08-16', '2021-08-17', '2021-08-18', '2021-08-19', '2021-08-20', '2021-08-21', '2021-08-22', '2021-08-23', '2021-08-24', '2021-08-25', '2021-08-26', '2021-08-27', '2021-08-28', '2021-08-29', '2021-08-30', '2021-08-31', '2021-09-01', '2021-09-02', '2021-09-03', '2021-09-04', '2021-09-05', '2021-09-06', '2021-09-07', '2021-09-08', '2021-09-09', '2021-09-10', '2021-09-11', '2021-09-12', '2021-09-13', '2021-09-14', '2021-09-15', '2021-09-16', '2021-09-17', '2021-09-18', '2021-09-19', '2021-09-20', '2021-09-21', '2021-09-22', '2021-09-23', '2021-09-24', '2021-09-25', '2021-09-26', '2021-09-27', '2021-09-28', '2021-09-29', '2021-09-30', '2021-10-01', '2021-10-02', '2021-10-03', '2021-10-04', '2021-10-05', '2021-10-06', '2021-10-07', '2021-10-08', '2021-10-09', '2021-10-10', '2021-10-11', '2021-10-12', '2021-10-13', '2021-10-14', '2021-10-15', '2021-10-16', '2021-10-17', '2021-10-18', '2021-10-19', '2021-10-20', '2021-10-21', '2021-10-22', '2021-10-23', '2021-10-24', '2021-10-25', '2021-10-26', '2021-10-27', '2021-10-28', '2021-10-29', '2021-10-30', '2021-10-31']
fig = plt.figure()
fig.set_size_inches(8,5.5)
ax = plt.subplot(1,1,1)
ax.set_xticks(linspace(0,len(dateIndex)-1,16))
ax.set_xticklabels(['01-21-20','03-04-20','04-16-20','05-29-20','07-12-20','08-24-20','10-06-20','11-18-20','01-01-21','02-13-21','03-28-21','05-10-21','06-23-21','08-05-21','09-17-21','10-31-21'],rotation=90)
ax.set_xlim((0,len(dateIndex)-1))
ax.fill_between(tSpan, susceptible_and_unvaccinated, facecolor = 'blue', label = 'Susceptible and unvaccinated', zorder = 5)
ax.fill_between(tSpan, susceptible_and_vaccinated + susceptible_and_unvaccinated, susceptible_and_unvaccinated, facecolor = 'green', label = 'Susceptible and vaccinated', zorder = 4)
ax.fill_between(tSpan, actively_infected + susceptible_and_vaccinated + susceptible_and_unvaccinated, susceptible_and_vaccinated + susceptible_and_unvaccinated, facecolor = 'orange', label = 'Actively infected', zorder = 3)
ax.fill_between(tSpan, actively_infected + susceptible_and_vaccinated + susceptible_and_unvaccinated + removed_and_unvaccinated, actively_infected + susceptible_and_vaccinated + susceptible_and_unvaccinated, facecolor = 'gray', label = 'Removed and unvaccinated', zorder = 2)
ax.fill_between(tSpan, removed_and_vaccinated + actively_infected + susceptible_and_vaccinated + susceptible_and_unvaccinated + removed_and_unvaccinated, actively_infected + susceptible_and_vaccinated + susceptible_and_unvaccinated + removed_and_unvaccinated, facecolor = 'yellow', label = 'Removed and vaccinated', zorder = 1)
ax.scatter(tSpan[::mult], cumsum(vr_daily[locationIndex]*S0), 10, marker = 'o', color = 'purple', label = 'Cumulative vaccinated', zorder = 6)
ax.set_ylabel('Proportion of initially susceptible population')
ax.set_ylim((0,S0))
ax.set_yticks(linspace(0,S0,5))
ax.set_yticklabels(['0','0.25','0.50','0.75','1'])
fig.tight_layout()
fig.legend(loc = (0.25,0.6), frameon=True, fontsize=8)
plt.savefig('Fig5C.pdf')
plt.close()