import numpy as np
import bionetgen
parLog = np.zeros((10000,20,18))
for ii in range(20):
    parLog[:,ii,:] = np.genfromtxt('/Users/amallela/Documents/perelson/sampling/output/Results/A_MCMC/Runs/params_'+str(ii)+'.txt')[33000:43000,:]    
parLog = parLog.reshape((200000,18))
scoreLog = np.genfromtxt('log.txt')[1160000:1360000,:]
scoreLog = scoreLog[np.lexsort(scoreLog.T[::-1])]
scoreLog = scoreLog[:,-1]
model = bionetgen.bngmodel("test.bngl")
kk = 1000
N = 58
p1,p2,p3,p4,p5,p6 = np.zeros((kk,561)),np.zeros((kk,561)),np.zeros((kk,561)),np.zeros((kk,561)),np.zeros((kk,561)),np.zeros((kk,561))
p7,p8,p9,p10 = np.zeros((kk,881)),np.zeros((kk,881)),np.zeros((kk,881)),np.zeros((kk,881))
sd = 0.15369023929
count = -1

for ii in np.arange(0,200000,200000//kk):
    count += 1
#     model.parameters.KaA0__FREE = parLog[ii,0]
#     model.parameters.KxRT__FREE = parLog[ii,1]
#     model.parameters.R1__FREE = parLog[ii,2]
#     model.parameters.R2__FREE = parLog[ii,3]
#     model.parameters.c__FREE = parLog[ii,4]
#     model.parameters.f__FREE = parLog[ii,5]
#     model.parameters.fp__FREE = parLog[ii,6]
#     model.parameters.fr__FREE = parLog[ii,7]
#     model.parameters.hatNc_ref__FREE = parLog[ii,8]
#     model.parameters.hatRT_ref__FREE = parLog[ii,9]
#     model.parameters.hata_ref__FREE = parLog[ii,10]
#     model.parameters.hatmu__FREE = parLog[ii,11]
#     model.parameters.hatp__FREE = parLog[ii,12]
#     model.parameters.kE__FREE = parLog[ii,13]
#     model.parameters.kfRT__FREE = parLog[ii,14]
#     model.parameters.kmc__FREE = parLog[ii,15]
#     model.parameters.kpc__FREE = parLog[ii,16]
#     model.parameters.rho__FREE = parLog[ii,17]
# 
#     with open("test_{}.bngl".format(ii), "w") as f:
#         f.write(str(model))
# 
#     with open("test_{}.bngl".format(ii), "r") as f:
#         fdata = f.read()
#         fdata = fdata.replace('Concentrations([])','Concentrations()')
#     with open("test_{}.bngl".format(ii), "w") as f:
#         f.write(fdata)

    result = bionetgen.run("test_"+str(ii)+".bngl")

    for jj in np.arange(561):
        p1[count,jj] = result["test_"+str(ii)+"_p1_control"][jj][16]
        p2[count,jj] = result["test_"+str(ii)+"_p1_control"][jj][21]
        p3[count,jj] = result["test_"+str(ii)+"_p1_control"][jj][24]
        p4[count,jj] = result["test_"+str(ii)+"_p1_treated"][jj][16]
        p5[count,jj] = result["test_"+str(ii)+"_p1_treated"][jj][21]
        p6[count,jj] = result["test_"+str(ii)+"_p1_treated"][jj][24]

    for jj in np.arange(881):
        p7[count,jj] = result["test_"+str(ii)+"_p2_control"][jj][16]
        p8[count,jj] = result["test_"+str(ii)+"_p2_control"][jj][29]
        p9[count,jj] = result["test_"+str(ii)+"_p2_treated"][jj][16]
        p10[count,jj] = result["test_"+str(ii)+"_p2_treated"][jj][29]

count = -1
for ii in np.arange(0,200000,200000//kk):
    count += 1
    p1[count,:] += np.random.normal(0.0, sd, 561)
    p2[count,:] += np.random.normal(0.0, sd, 561)
    p3[count,:] += np.random.normal(0.0, sd, 561)
    p4[count,:] += np.random.normal(0.0, sd, 561)
    p5[count,:] += np.random.normal(0.0, sd, 561)
    p6[count,:] += np.random.normal(0.0, sd, 561)
    p7[count,:] += np.random.normal(0.0, sd, 881)
    p8[count,:] += np.random.normal(0.0, sd, 881)
    p9[count,:] += np.random.normal(0.0, sd, 881)
    p10[count,:] += np.random.normal(0.0, sd, 881) 

p1l68 = np.percentile(p1, 16, axis=0)
p1u68 = np.percentile(p1, 84, axis=0)
p1l95 = np.percentile(p1, 2.5, axis=0)
p1u95 = np.percentile(p1, 97.5, axis=0)
p2l68 = np.percentile(p2, 16, axis=0)
p2u68 = np.percentile(p2, 84, axis=0)
p2l95 = np.percentile(p2, 2.5, axis=0)
p2u95 = np.percentile(p2, 97.5, axis=0)
p3l68 = np.percentile(p3, 16, axis=0)
p3u68 = np.percentile(p3, 84, axis=0)
p3l95 = np.percentile(p3, 2.5, axis=0)
p3u95 = np.percentile(p3, 97.5, axis=0)
p4l68 = np.percentile(p4, 16, axis=0)
p4u68 = np.percentile(p4, 84, axis=0)
p4l95 = np.percentile(p4, 2.5, axis=0)
p4u95 = np.percentile(p4, 97.5, axis=0)
p5l68 = np.percentile(p5, 16, axis=0)
p5u68 = np.percentile(p5, 84, axis=0)
p5l95 = np.percentile(p5, 2.5, axis=0)
p5u95 = np.percentile(p5, 97.5, axis=0)
p6l68 = np.percentile(p6, 16, axis=0)
p6u68 = np.percentile(p6, 84, axis=0)
p6l95 = np.percentile(p6, 2.5, axis=0)
p6u95 = np.percentile(p6, 97.5, axis=0)
p7l68 = np.percentile(p7, 16, axis=0)
p7u68 = np.percentile(p7, 84, axis=0)
p7l95 = np.percentile(p7, 2.5, axis=0)
p7u95 = np.percentile(p7, 97.5, axis=0)
p8l68 = np.percentile(p8, 16, axis=0)
p8u68 = np.percentile(p8, 84, axis=0)
p8l95 = np.percentile(p8, 2.5, axis=0)
p8u95 = np.percentile(p8, 97.5, axis=0)
p9l68 = np.percentile(p9, 16, axis=0)
p9u68 = np.percentile(p9, 84, axis=0)
p9l95 = np.percentile(p9, 2.5, axis=0)
p9u95 = np.percentile(p9, 97.5, axis=0)
p10l68 = np.percentile(p10, 16, axis=0)
p10u68 = np.percentile(p10, 84, axis=0)
p10l95 = np.percentile(p10, 2.5, axis=0)
p10u95 = np.percentile(p10, 97.5, axis=0)
p1_median = np.median(p1, axis=0)
p2_median = np.median(p2, axis=0)
p3_median = np.median(p3, axis=0)
p4_median = np.median(p4, axis=0)
p5_median = np.median(p5, axis=0)
p6_median = np.median(p6, axis=0)
p7_median = np.median(p7, axis=0)
p8_median = np.median(p8, axis=0)
p9_median = np.median(p9, axis=0)
p10_median = np.median(p10, axis=0)

np.save('p1l68.npy',p1l68)
np.save('p1l95.npy',p1l95)
np.save('p1u68.npy',p1u68)
np.save('p1u95.npy',p1u95)
np.save('p2l68.npy',p2l68)
np.save('p2l95.npy',p2l95)
np.save('p2u68.npy',p2u68)
np.save('p2u95.npy',p2u95)
np.save('p3l68.npy',p3l68)
np.save('p3l95.npy',p3l95)
np.save('p3u68.npy',p3u68)
np.save('p3u95.npy',p3u95)
np.save('p4l68.npy',p4l68)
np.save('p4l95.npy',p4l95)
np.save('p4u68.npy',p4u68)
np.save('p4u95.npy',p4u95)
np.save('p5l68.npy',p5l68)
np.save('p5l95.npy',p5l95)
np.save('p5u68.npy',p5u68)
np.save('p5u95.npy',p5u95)
np.save('p6l68.npy',p6l68)
np.save('p6l95.npy',p6l95)
np.save('p6u68.npy',p6u68)
np.save('p6u95.npy',p6u95)
np.save('p7l68.npy',p7l68)
np.save('p7l95.npy',p7l95)
np.save('p7u68.npy',p7u68)
np.save('p7u95.npy',p7u95)
np.save('p8l68.npy',p8l68)
np.save('p8l95.npy',p8l95)
np.save('p8u68.npy',p8u68)
np.save('p8u95.npy',p8u95)
np.save('p9l68.npy',p9l68)
np.save('p9l95.npy',p9l95)
np.save('p9u68.npy',p9u68)
np.save('p9u95.npy',p9u95)
np.save('p10l68.npy',p10l68)
np.save('p10l95.npy',p10l95)
np.save('p10u68.npy',p10u68)
np.save('p10u95.npy',p10u95)
np.save('p1_median.npy',p1_median)
np.save('p2_median.npy',p2_median)
np.save('p3_median.npy',p3_median)
np.save('p4_median.npy',p4_median)
np.save('p5_median.npy',p5_median)
np.save('p6_median.npy',p6_median)
np.save('p7_median.npy',p7_median)
np.save('p8_median.npy',p8_median)
np.save('p9_median.npy',p9_median)
np.save('p10_median.npy',p10_median)