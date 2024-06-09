import numpy as np
import pickle
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':3})
with open('/Users/amallela/Documents/codes/EID/ADJ/LSODA/withJacobian/adjy','rb') as f1:
    adjy = pickle.load(f1)
with open('/Users/amallela/Documents/codes/EID/CFD/LSODA/withJacobian/cfdy','rb') as f3:
    cfdy = pickle.load(f3)        
ady = np.load('/Users/amallela/Documents/codes/EID/AD/DOPRI8/ss.npy')
fwdy = np.array([ 1.21384702e+00,  1.86644932e+01,  8.56555970e-01,  1.62656599e-02,  1.14809785e+02, -2.88624083e+03, -1.45783939e+03,  9.40671724e-02,  7.93015300e+02, -3.56378069e-03, -1.38152284e+01, -1.14756312e-05,  3.85283177e+02,  1.83925658e+02, -3.88212882e+00,  1.04469777e+02, -3.57252909e+00,  2.43835599e+01,  0.00000000e+00, -2.26039197e+01, -3.57252909e+00,  7.72779565e+01,  0.00000000e+00,  0.00000000e+00, -1.82713695e+02,  2.31191765e-01])
parLabels = np.array(['$t_0$', '$\\sigma - t_0$', '$\\tau_1 - \\sigma$', '$\\tau_2 - \\tau_1$','$\\beta$','$\\lambda_0$','$p_0$','$\\lambda_1$','$p_1$','$\\lambda_2$','$p_2$','$S_0$','$m_b$','$\\rho_{\\rm E}$','$\\rho_{\\rm A}$','$k_L$','$c_I$','$f_A$','$f_H$','$k_Q$','$j_Q$','$c_A$','$c_H$','$f_R$','$f_D$','$r$'])
plt.figure(figsize=(100,100))
fig, ax = plt.subplots(nrows=5,ncols=5)
count = -1
for i in range(5):
    for j in range(5):
        count += 1
        if count < 22:
            ax[i,j].scatter(np.arange(1), ady[i*5+j+4][0], 4, marker='s', facecolor=[], edgecolor='k', label='AD, 1e-8')
            ax[i,j].scatter(np.arange(1), ady[i*5+j+4][2], 4, marker='o', facecolor=[], edgecolor='y', label='AD, 1e-6')            
            ax[i,j].scatter(np.arange(1), cfdy[-1][0][i*5+j+4], 4, marker='^', facecolor=[], edgecolor='g', label='FD')
            ax[i,j].scatter(np.arange(1), adjy[-1][0][i*5+j+4], 4, marker='v', facecolor=[], edgecolor='r', label='ADJ')
            ax[i,j].scatter(np.arange(1), fwdy[i*5+j+4], 4, marker='d', facecolor=[], edgecolor='b', label='FWD')
            ax[i,j].set_xticklabels([])
            ax[i,j].set_xlabel(parLabels[i*5+j+4],fontsize=4)
            ax[i,j].ticklabel_format(axis='y', style='sci')
ax[2,2].legend(frameon=False)
plt.tight_layout(pad=5.0)
plt.savefig('parameter_plots.pdf')