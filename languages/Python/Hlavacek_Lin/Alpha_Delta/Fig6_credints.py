import numpy as np
import arviz as az

arr_alpha_Dallas, arr_delta_Dallas, arr_alpha_Phoenix, arr_delta_Phoenix, arr_alpha_Houston, arr_delta_Houston, arr_alpha_NYC, arr_delta_NYC = np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3)

l = np.load('inferences/Dallas4-n5-logLLog.npy')
x = np.load('inferences/Dallas4-n5-parLog.npy')
bestFitPar = x[np.where(l==np.amax(l))[0][0]]
arr_alpha_Dallas[0] = bestFitPar[16]
arr_alpha_Dallas[1:3] = az.hdi(x[220000:340000,16],hdi_prob=0.95)
arr_delta_Dallas[0] = bestFitPar[16] + bestFitPar[17]
arr_delta_Dallas[1:3] = bestFitPar[16] + az.hdi(x[220000:340000,17],hdi_prob=0.95)

l = np.load('inferences/Phoenix29-n5-logLLog.npy')
x = np.load('inferences/Phoenix29-n5-parLog.npy')
bestFitPar = x[np.where(l==np.amax(l))[0][0]]
arr_alpha_Phoenix[0] = bestFitPar[16]
arr_alpha_Phoenix[1:3] = az.hdi(x[0:100000,16],hdi_prob=0.95)
arr_delta_Phoenix[0] = bestFitPar[16] + bestFitPar[17]
arr_delta_Phoenix[1:3] = bestFitPar[16] + az.hdi(x[0:100000,17],hdi_prob=0.95)

l = np.load('inferences/Houston2-n4-logLLog.npy')
x = np.load('inferences/Houston2-n4-parLog.npy')
bestFitPar = x[np.where(l==np.amax(l))[0][0]]
arr_alpha_Houston[0] = bestFitPar[13]
arr_alpha_Houston[1:3] = az.hdi(x[200000:380000,13],hdi_prob=0.95)
arr_delta_Houston[0] = bestFitPar[13] + bestFitPar[14]
arr_delta_Houston[1:3] = bestFitPar[13] + az.hdi(x[200000:380000,14],hdi_prob=0.95)

l = np.load('inferences/NYC0-n4-logLLog.npy')
x = np.load('inferences/NYC0-n4-parLog.npy')
bestFitPar = x[np.where(l==np.amax(l))[0][0]]
arr_alpha_NYC[0] = bestFitPar[13]
arr_alpha_NYC[1:3] = az.hdi(x[12000:36000,13],hdi_prob=0.95)
arr_delta_NYC[0] = bestFitPar[13] + bestFitPar[14]
arr_delta_NYC[1:3] = bestFitPar[13] + az.hdi(x[12000:36000,14],hdi_prob=0.95)

np.save('arr_alpha_Dallas.npy',arr_alpha_Dallas)
np.save('arr_delta_Dallas.npy',arr_delta_Dallas)
np.save('arr_alpha_Phoenix.npy',arr_alpha_Phoenix)
np.save('arr_delta_Phoenix.npy',arr_delta_Phoenix)
np.save('arr_alpha_Houston.npy',arr_alpha_Houston)
np.save('arr_delta_Houston.npy',arr_delta_Houston)
np.save('arr_alpha_NYC.npy',arr_alpha_NYC)
np.save('arr_delta_NYC.npy',arr_delta_NYC)