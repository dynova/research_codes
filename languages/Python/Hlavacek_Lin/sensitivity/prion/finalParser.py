from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from scipy.integrate import solve_ivp
from scipy.integrate import ode
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import time

m = open('SKMEL28ode.net','r')
#m = open('testModel.net','r')

mLines = m.readlines()

for i in range(len(mLines)):
    
    if ('begin' in mLines[i])&('parameters' in mLines[i]):
        parInit = i+1
        
    if ('end' in mLines[i])&('parameters' in mLines[i]):
        parEnd = i
    
    if ('begin' in mLines[i])&('reactions' in mLines[i]):
        reactInit = i+1
        
    if ('end' in mLines[i])&('reactions' in mLines[i]):
        reactEnd = i
        
    if ('begin' in mLines[i])&('species' in mLines[i]):
        speciesInit = i+1
        
    if ('end' in mLines[i])&('species' in mLines[i]):
        speciesEnd = i    
    
parLines = mLines[parInit:parEnd]
lines = mLines[reactInit:reactEnd]  
speciesLines = mLines[speciesInit:speciesEnd]  
numPars = len(parLines)

for parID in range(len(parLines)):
    
    parLines[parID]=parLines[parID].replace("^","**")
    
    separated= parLines[parID].split(' ')
    
    for index in range(len(separated)-1,-1,-1):
        
        if len(separated[index])==0:
            del separated[index]

    exec(separated[1]+'='+separated[2])

numSpecies = len(speciesLines)
IC= zeros((numSpecies,))

for speciesID in range(numSpecies):
    
    separated= speciesLines[speciesID].split(' ')
    
    for index in range(len(separated)-1,-1,-1):
        
        if len(separated[index])==0:
            del separated[index]

    exec('IC[speciesID]='+separated[2])


reactants = zeros((len(lines),10))
products = zeros((len(lines),10))

rates = zeros((len(lines), ))

for reactionID in range(len(lines)):
    
    separated= lines[reactionID].split(' ')
    
    for index in range(len(separated)-1,-1,-1):
        
        if len(separated[index])==0:
            del separated[index]
            
    reactantSet = separated[1].split(',')
    reactants[reactionID][0]=len(reactantSet)
    
    for reactantID in range(len(reactantSet)):
        reactants[reactionID][reactantID+1] = int(reactantSet[reactantID])-1
        
    productSet = separated[2].split(',')
    products[reactionID][0]=len(productSet)
    
    for productID in range(len(productSet)):
        products[reactionID][productID+1] = int(productSet[productID])-1
    
#    rates[reactionID].append(eval(separated[3])) 
    rates[reactionID]=eval(separated[3])


numReactions = len(rates)

minNum = 20
dimerizationIndex = zeros((numReactions,))
rdimerizationIndex = zeros((numReactions,))

for i in range(len(rates)):

    if minNum > reactants[i][0]:
        
        minNum = reactants[i][0]
    
    if (reactants[i][0]==2)&(reactants[i][1]==reactants[i][2]):
        print('reaction ID='+str(i)+' is dimerization!!')
        print(asarray(reactants[i])+1)
        print(asarray(products[i])+1)
        
        dimerizationIndex[i]=1.

reactants=reactants.astype('int')
products=products.astype('int')
        
def dX(t,X):

    dx = zeros((numSpecies,))
    
    for i in range(numReactions):

        flux = rates[i]
        
        for j in range(1, reactants[i][0]+1):
            
            flux = flux*X[reactants[i][j]]
        
        for j in range(1, reactants[i][0]+1):

            dx[ reactants[i][j] ]  -= flux
            
        for j in range(1, products[i][0]+1):

            dx[ products[i][j] ]  += flux

    return dx


tStart = time.time()
method='LSODA'
tSpan = linspace(0,5,1000)
sol = solve_ivp(dX,(0,tSpan[-1]),IC,t_eval=tSpan,method=method,rtol=1e-3,atol=1e-6)
tEnd = time.time()

print(tEnd-tStart)

indexForPlotting = [186,165,6]
#indexForPlotting = [0,1]

fig = plt.figure()
fig.set_size_inches((15,7))

checks = genfromtxt('generic.cdat')

for i in indexForPlotting:

    plt.plot(sol.t, sol.y[i-1,:], lw=2,color='b')    
    plt.scatter(checks[:,0], checks[:,i], s=50,color='r')    

#plt.scatter(range(numSpecies),sol.y[:,0], s=50, color='g')
#plt.scatter(range(numSpecies),checks[0,1:], s=20, color='r')