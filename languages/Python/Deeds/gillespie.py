############################################################
#### Gillespie algorithm for ProcProc Model                #
#### Written by Abhishek Mallela                           #
#### Pseudocode + algorithm developed with Eric J. Deeds   #
#### Latest Version: 8/18/16                               #
############################################################

import numpy as np
from numpy import *
import random
import scipy.io

Mon = 1e-5 # This is for r = ratio of V_max's = 1e-2
Don = 1e-5 # Stochastic 'on' rate
r1 = 2e-5
r2 = 2e-4
Moff = 1e-1 # This is for r = ratio of V_max's = 1e-2; keeping K_M the same
Doff = 1e-1 # Stochastic, deterministic: same
Mcat2 = 1-1e-3 # This is for r = ratio of V_max's = 1e-2; keeping K_M the same
Mcat1 = Mcat2/100 # Because it takes longer to add the first ubiquitin unit
Dcat = 1-1e-3 # Stochastic, deterministic: same

def gillespie(maxtime=1e6):
    #### Initial declarations
    num = 1000 # S_init
    c = 0
    Q = r1*num # Defining it this way so that 'num' is the only quantity to change
    v0 = np.zeros((num)) # Row vector of 'num'-many zeros. Each element represents an S0 molecule
    v1 = np.empty((c)) # Empty row vector for v1; container for S1 to S3
    v2 = np.empty((c)) # Empty row vector for v2; container for S4 to Sl
    v3 = np.empty((c)) # Empty row vector for v3; container for MS0 to MS3
    v4 = np.empty((c)) # Empty row vector for v4; container for DS0 to DS3
    v5 = np.empty((c)) # Empty row vector for v5; container for MS4 to MSl
    v6 = np.empty((c)) # Empty row vector for v6; container for DS4 to DSl
    M_un = 100 # M_init
    D_un = 100 # D_init
    
## Begin Stochastic Simulation Algorithm (Doob-Gillespie)
    time = 0
    while (time < maxtime):
        tot = [Q, r1*len(v1), r2*len(v2), r1*len(v3), r1*len(v4), + 
    					r2*len(v5), r2*len(v6), Mon*M_un*len(v1), +
    					Don*D_un*len(v1), Mon*M_un*len(v0), Don*D_un*len(v0), + 
    					Mon*M_un*len(v2), Don*D_un*len(v2), Moff*len(v3), + 
    					Doff*len(v4), Moff*len(v5), Doff*len(v6), +
    					Mcat1*(len(v3)-np.count_nonzero(v3)), +
    					Mcat2*np.count_nonzero(v3), +
    					Dcat*len(v4), Mcat2*len(v5), Dcat*len(v6)]
        for k in range(21):
            tot[k + 1] += tot[k]
        kT = tot[-1]
        time += -np.log(1-np.random.random())/kT
        
        if time >= maxtime:
            break        
            
        event_kT = np.random.random() * kT
        #### Reaction 1: Synthesize a new S
        if event_kT < tot[0]:
            v0 = np.append(v0,0)
        #### Reaction 2: Delete a S from v1
        elif event_kT < tot[1]:
            k = random.randrange(len(v1))
            v1 = v1[np.arange(len(v1))!=k]
        #### Reaction 3: Delete a S from v2        
        elif event_kT < tot[2]:
            k = random.randrange(len(v2))
            v2 = v2[np.arange(len(v2))!=k]
        #### Reaction 4: Delete a S from v3 
        elif event_kT < tot[3]:
            k = random.randrange(len(v3))
            v3 = v3[np.arange(len(v3))!=k]
            M_un += 1
        #### Reaction 5: Delete a S from v4 
        elif event_kT < tot[4]:
            k = random.randrange(len(v4))
            v4 = v4[np.arange(len(v4))!=k]
            D_un += 1
        #### Reaction 6: Delete a S from v5
        elif event_kT < tot[5]:
            k = random.randrange(len(v5))
            v5 = v5[np.arange(len(v5))!=k]
            M_un += 1
        #### Reaction 7: Delete a S from v6
        elif event_kT < tot[6]:
            k = random.randrange(len(v6))
            v6 = v6[np.arange(len(v6))!=k]
            D_un += 1
        #### Reaction 8: Bind a M to an S in v1
        elif event_kT < tot[7]:
            k = random.randrange(len(v1))
            v1temp = v1[np.arange(len(v1))!=k]
            v1rem = v1[k]
            v1 = v1temp
            if (v1rem > 0 and v1rem < 4):
                v3 = np.append(v3,v1rem)
            elif (v1rem >= 4):
                v5 = np.append(v5,v1rem)
            M_un -= 1
        #### Reaction 9: Bind a D to an S in v1
        elif event_kT < tot[8]:
            k = random.randrange(len(v1))
            v1temp = v1[np.arange(len(v1))!=k]
            v1rem = v1[k]
            v1 = v1temp
            if (v1rem > 0 and v1rem < 4):
                v4 = np.append(v4,v1rem)
            elif (v1rem >= 4):
                v6 = np.append(v6,v1rem)
            D_un -= 1
        #### Reaction 10: Bind a M to an S in v0
        elif event_kT < tot[9]:
            k = random.randrange(len(v0))
            v0temp = v0[np.arange(len(v0))!=k]
            v0rem = v0[k]
            v0 = v0temp
            v3 = np.append(v3,v0rem)
            M_un -= 1
        #### Reaction 11: Bind a D to an S in v0
        elif event_kT < tot[10]:
            k = random.randrange(len(v0))
            v0temp = v0[np.arange(len(v0))!=k]
            v0rem = v0[k]
            v0 = v0temp
            v4 = np.append(v4,v0rem)
            D_un -= 1
        #### Reaction 12: Bind a M to an S in v2
        elif event_kT < tot[11]:
            k = random.randrange(len(v2))
            v2temp = v2[np.arange(len(v2))!=k]
            v2rem = v2[k]
            v2 = v2temp
            v5 = np.append(v5,v2rem)
            M_un -= 1
        #### Reaction 13: Bind a D to an S in v2
        elif event_kT < tot[12]:
            k = random.randrange(len(v2))
            v2temp = v2[np.arange(len(v2))!=k]
            v2rem = v2[k]
            v2 = v2temp
            v6 = np.append(v6,v2rem)
            D_un -= 1
        #### Reaction 14: Unbind a M from v3
        elif event_kT < tot[13]:
            k = random.randrange(len(v3))
            v3temp = v3[np.arange(len(v3))!=k]
            v3rem = v3[k]
            if v3rem == 0:
                v0 = np.append(v0,v3rem)
            elif v3rem > 0:
                v1 = np.append(v1,v3rem)
            v3 = v3temp
            M_un += 1
        #### Reaction 15: Unbind a D from v4
        elif event_kT < tot[14]:
            k = random.randrange(len(v4))
            v4temp = v4[np.arange(len(v4))!=k]
            v4rem = v4[k]
            if v4rem == 0:
                v0 = np.append(v0,v4rem)
            elif v4rem > 0:
                v1 = np.append(v1,v4rem)
            v4 = v4temp
            D_un += 1
        #### Reaction 16: Unbind a M from v5
        elif event_kT < tot[15]:
            k = random.randrange(len(v5))
            v5temp = v5[np.arange(len(v5))!=k]
            v5rem = v5[k]
            v2 = np.append(v2,v5rem)
            v5 = v5temp
            M_un += 1
        #### Reaction 17: Unbind a D from v6
        elif event_kT < tot[16]:
            k = random.randrange(len(v6))
            v6temp = v6[np.arange(len(v6))!=k]
            v6rem = v6[k]
            v2 = np.append(v2,v6rem)
            v6 = v6temp
            D_un += 1
        #### Reaction 18: M catalysis for Mcat1 reaction in v3
        elif event_kT < tot[17]:
            k = random.randrange(len(v3))
            if v3[k] == 0:
                v3[k] += 1
        #### Reaction 19: M catalysis for Mcat2 reactions in v3   
        elif event_kT < tot[18]:
            k = random.randrange(len(v3))
            if v3[k] == 1 or v3[k] == 2 or v3[k] == 3:
                v3[k] += 1
                if v3[k] > 3:
                    v5 = np.append(v5,v3[k])
                    v3 = v3[np.arange(len(v3))!=k]
        #### Reaction 20: D catalysis for v4
        elif event_kT < tot[19]:
            k = random.randrange(len(v4))
            v4[k] -= 1
            if v4[k] == 0:
                v4 = v4[np.arange(len(v4))!=k]
                v0 = np.append(v0,0)
        #### Reaction 21: M catalysis for v5
        elif event_kT < tot[20]:
            k = random.randrange(len(v5))
            v5[k] += 1
        #### Reaction 22: D catalysis for v6
        elif event_kT < tot[21]:
            k = random.randrange(len(v6))
            v6[k] -= 1
            if v6[k] == 3:
                v4 = np.append(v4,v6[k])        
                v6 = v6[np.arange(len(v6))!=k]
     
    l0 = len(v0)
    l1 = len(v1)
    l2 = len(v2)
    l3 = len(v3)
    l4 = len(v4)
    l5 = len(v5)
    l6 = len(v6)
    l7 = max(v2)
    return l0, l1, l2, l3, l4, l5, l6, l7

if __name__ == '__main__':
    a = gillespie()
    scipy.io.savemat('ProcProc1.mat',dict(Smod=a[2]/10,Stot=sum(a)/10 - a[7]/10,lmax = a[7]))