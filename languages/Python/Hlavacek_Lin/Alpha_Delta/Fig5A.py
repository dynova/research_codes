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
from scipy.io import loadmat
plt.rcParams.update({'font.family':'sans','font.size':12,'image.cmap':'plasma'})
colors1 = cm.Blues(linspace(0,1,9))
colors2 = cm.Reds(linspace(0,1,9))
dateIndex = ['2020-01-21', '2020-01-22', '2020-01-23', '2020-01-24', '2020-01-25', '2020-01-26', '2020-01-27', '2020-01-28', '2020-01-29', '2020-01-30', '2020-01-31', '2020-02-01', '2020-02-02', '2020-02-03', '2020-02-04', '2020-02-05', '2020-02-06', '2020-02-07', '2020-02-08', '2020-02-09', '2020-02-10', '2020-02-11', '2020-02-12', '2020-02-13', '2020-02-14', '2020-02-15', '2020-02-16', '2020-02-17', '2020-02-18', '2020-02-19', '2020-02-20', '2020-02-21', '2020-02-22', '2020-02-23', '2020-02-24', '2020-02-25', '2020-02-26', '2020-02-27', '2020-02-28', '2020-02-29', '2020-03-01', '2020-03-02', '2020-03-03', '2020-03-04', '2020-03-05', '2020-03-06', '2020-03-07', '2020-03-08', '2020-03-09', '2020-03-10', '2020-03-11', '2020-03-12', '2020-03-13', '2020-03-14', '2020-03-15', '2020-03-16', '2020-03-17', '2020-03-18', '2020-03-19', '2020-03-20', '2020-03-21', '2020-03-22', '2020-03-23', '2020-03-24', '2020-03-25', '2020-03-26', '2020-03-27', '2020-03-28', '2020-03-29', '2020-03-30', '2020-03-31', '2020-04-01', '2020-04-02', '2020-04-03', '2020-04-04', '2020-04-05', '2020-04-06', '2020-04-07', '2020-04-08', '2020-04-09', '2020-04-10', '2020-04-11', '2020-04-12', '2020-04-13', '2020-04-14', '2020-04-15', '2020-04-16', '2020-04-17', '2020-04-18', '2020-04-19', '2020-04-20', '2020-04-21', '2020-04-22', '2020-04-23', '2020-04-24', '2020-04-25', '2020-04-26', '2020-04-27', '2020-04-28', '2020-04-29', '2020-04-30', '2020-05-01', '2020-05-02', '2020-05-03', '2020-05-04', '2020-05-05', '2020-05-06', '2020-05-07', '2020-05-08', '2020-05-09', '2020-05-10', '2020-05-11', '2020-05-12', '2020-05-13', '2020-05-14', '2020-05-15', '2020-05-16', '2020-05-17', '2020-05-18', '2020-05-19', '2020-05-20', '2020-05-21', '2020-05-22', '2020-05-23', '2020-05-24', '2020-05-25', '2020-05-26', '2020-05-27', '2020-05-28', '2020-05-29', '2020-05-30', '2020-05-31', '2020-06-01', '2020-06-02', '2020-06-03', '2020-06-04', '2020-06-05', '2020-06-06', '2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', '2020-06-12', '2020-06-13', '2020-06-14', '2020-06-15', '2020-06-16', '2020-06-17', '2020-06-18', '2020-06-19', '2020-06-20', '2020-06-21', '2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25', '2020-06-26', '2020-06-27', '2020-06-28', '2020-06-29', '2020-06-30', '2020-07-01', '2020-07-02', '2020-07-03', '2020-07-04', '2020-07-05', '2020-07-06', '2020-07-07', '2020-07-08', '2020-07-09', '2020-07-10', '2020-07-11', '2020-07-12', '2020-07-13', '2020-07-14', '2020-07-15', '2020-07-16', '2020-07-17', '2020-07-18', '2020-07-19', '2020-07-20', '2020-07-21', '2020-07-22', '2020-07-23', '2020-07-24', '2020-07-25', '2020-07-26', '2020-07-27', '2020-07-28', '2020-07-29', '2020-07-30', '2020-07-31', '2020-08-01', '2020-08-02', '2020-08-03', '2020-08-04', '2020-08-05', '2020-08-06', '2020-08-07', '2020-08-08', '2020-08-09', '2020-08-10', '2020-08-11', '2020-08-12', '2020-08-13', '2020-08-14', '2020-08-15', '2020-08-16', '2020-08-17', '2020-08-18', '2020-08-19', '2020-08-20', '2020-08-21', '2020-08-22', '2020-08-23', '2020-08-24', '2020-08-25', '2020-08-26', '2020-08-27', '2020-08-28', '2020-08-29', '2020-08-30', '2020-08-31', '2020-09-01', '2020-09-02', '2020-09-03', '2020-09-04', '2020-09-05', '2020-09-06', '2020-09-07', '2020-09-08', '2020-09-09', '2020-09-10', '2020-09-11', '2020-09-12', '2020-09-13', '2020-09-14', '2020-09-15', '2020-09-16', '2020-09-17', '2020-09-18', '2020-09-19', '2020-09-20', '2020-09-21', '2020-09-22', '2020-09-23', '2020-09-24', '2020-09-25', '2020-09-26', '2020-09-27', '2020-09-28', '2020-09-29', '2020-09-30', '2020-10-01', '2020-10-02', '2020-10-03', '2020-10-04', '2020-10-05', '2020-10-06', '2020-10-07', '2020-10-08', '2020-10-09', '2020-10-10', '2020-10-11', '2020-10-12', '2020-10-13', '2020-10-14', '2020-10-15', '2020-10-16', '2020-10-17', '2020-10-18', '2020-10-19', '2020-10-20', '2020-10-21', '2020-10-22', '2020-10-23', '2020-10-24', '2020-10-25', '2020-10-26', '2020-10-27', '2020-10-28', '2020-10-29', '2020-10-30', '2020-10-31', '2020-11-01', '2020-11-02', '2020-11-03', '2020-11-04', '2020-11-05', '2020-11-06', '2020-11-07', '2020-11-08', '2020-11-09', '2020-11-10', '2020-11-11', '2020-11-12', '2020-11-13', '2020-11-14', '2020-11-15', '2020-11-16', '2020-11-17', '2020-11-18', '2020-11-19', '2020-11-20', '2020-11-21', '2020-11-22', '2020-11-23', '2020-11-24', '2020-11-25', '2020-11-26', '2020-11-27', '2020-11-28', '2020-11-29', '2020-11-30','2020-12-01', '2020-12-02', '2020-12-03', '2020-12-04', '2020-12-05', '2020-12-06', '2020-12-07', '2020-12-08', '2020-12-09', '2020-12-10', '2020-12-11', '2020-12-12', '2020-12-13', '2020-12-14', '2020-12-15', '2020-12-16', '2020-12-17', '2020-12-18', '2020-12-19', '2020-12-20', '2020-12-21', '2020-12-22', '2020-12-23', '2020-12-24', '2020-12-25', '2020-12-26', '2020-12-27', '2020-12-28', '2020-12-29', '2020-12-30','2020-12-31', '2021-01-01', '2021-01-02', '2021-01-03', '2021-01-04', '2021-01-05', '2021-01-06', '2021-01-07', '2021-01-08', '2021-01-09', '2021-01-10', '2021-01-11', '2021-01-12', '2021-01-13', '2021-01-14', '2021-01-15', '2021-01-16', '2021-01-17', '2021-01-18', '2021-01-19', '2021-01-20', '2021-01-21', '2021-01-22', '2021-01-23', '2021-01-24', '2021-01-25', '2021-01-26', '2021-01-27', '2021-01-28', '2021-01-29', '2021-01-30', '2021-01-31', '2021-02-01', '2021-02-02', '2021-02-03', '2021-02-04', '2021-02-05', '2021-02-06', '2021-02-07', '2021-02-08', '2021-02-09', '2021-02-10', '2021-02-11', '2021-02-12', '2021-02-13', '2021-02-14', '2021-02-15', '2021-02-16', '2021-02-17', '2021-02-18', '2021-02-19', '2021-02-20', '2021-02-21', '2021-02-22', '2021-02-23', '2021-02-24', '2021-02-25', '2021-02-26', '2021-02-27', '2021-02-28', '2021-02-29', '2021-03-01', '2021-03-02', '2021-03-03', '2021-03-04', '2021-03-05', '2021-03-06', '2021-03-07', '2021-03-08', '2021-03-09', '2021-03-10', '2021-03-11', '2021-03-12', '2021-03-13', '2021-03-14', '2021-03-15', '2021-03-16', '2021-03-17', '2021-03-18', '2021-03-19', '2021-03-20', '2021-03-21', '2021-03-22', '2021-03-23', '2021-03-24', '2021-03-25', '2021-03-26', '2021-03-27', '2021-03-28', '2021-03-29', '2021-03-30', '2021-03-31', '2021-04-01', '2021-04-02', '2021-04-03', '2021-04-04', '2021-04-05', '2021-04-06', '2021-04-07', '2021-04-08', '2021-04-09', '2021-04-10', '2021-04-11', '2021-04-12', '2021-04-13', '2021-04-14', '2021-04-15', '2021-04-16', '2021-04-17', '2021-04-18', '2021-04-19', '2021-04-20', '2021-04-21', '2021-04-22', '2021-04-23', '2021-04-24', '2021-04-25', '2021-04-26', '2021-04-27', '2021-04-28', '2021-04-29', '2021-04-30', '2021-05-01', '2021-05-02', '2021-05-03', '2021-05-04', '2021-05-05', '2021-05-06', '2021-05-07', '2021-05-08', '2021-05-09', '2021-05-10', '2021-05-11', '2021-05-12', '2021-05-13', '2021-05-14', '2021-05-15', '2021-05-16', '2021-05-17', '2021-05-18', '2021-05-19', '2021-05-20', '2021-05-21', '2021-05-22', '2021-05-23', '2021-05-24', '2021-05-25', '2021-05-26', '2021-05-27', '2021-05-28', '2021-05-29', '2021-05-30', '2021-05-31', '2021-06-01', '2021-06-02', '2021-06-03', '2021-06-04', '2021-06-05', '2021-06-06', '2021-06-07', '2021-06-08', '2021-06-09', '2021-06-10', '2021-06-11', '2021-06-12', '2021-06-13', '2021-06-14', '2021-06-15', '2021-06-16', '2021-06-17', '2021-06-18', '2021-06-19', '2021-06-20', '2021-06-21', '2021-06-22', '2021-06-23', '2021-06-24', '2021-06-25', '2021-06-26', '2021-06-27', '2021-06-28', '2021-06-29', '2021-06-30', '2021-07-01', '2021-07-02', '2021-07-03', '2021-07-04', '2021-07-05', '2021-07-06', '2021-07-07', '2021-07-08', '2021-07-09', '2021-07-10', '2021-07-11', '2021-07-12', '2021-07-13', '2021-07-14', '2021-07-15', '2021-07-16', '2021-07-17', '2021-07-18', '2021-07-19', '2021-07-20', '2021-07-21', '2021-07-22', '2021-07-23', '2021-07-24', '2021-07-25', '2021-07-26', '2021-07-27', '2021-07-28', '2021-07-29', '2021-07-30', '2021-07-31', '2021-08-01', '2021-08-02', '2021-08-03', '2021-08-04', '2021-08-05', '2021-08-06', '2021-08-07', '2021-08-08', '2021-08-09', '2021-08-10', '2021-08-11', '2021-08-12', '2021-08-13', '2021-08-14', '2021-08-15', '2021-08-16', '2021-08-17', '2021-08-18', '2021-08-19', '2021-08-20', '2021-08-21', '2021-08-22', '2021-08-23', '2021-08-24', '2021-08-25', '2021-08-26', '2021-08-27', '2021-08-28', '2021-08-29', '2021-08-30', '2021-08-31', '2021-09-01', '2021-09-02', '2021-09-03', '2021-09-04', '2021-09-05', '2021-09-06', '2021-09-07', '2021-09-08', '2021-09-09', '2021-09-10', '2021-09-11', '2021-09-12', '2021-09-13', '2021-09-14', '2021-09-15', '2021-09-16', '2021-09-17', '2021-09-18', '2021-09-19', '2021-09-20', '2021-09-21', '2021-09-22', '2021-09-23', '2021-09-24', '2021-09-25', '2021-09-26', '2021-09-27', '2021-09-28', '2021-09-29', '2021-09-30', '2021-10-01', '2021-10-02', '2021-10-03', '2021-10-04', '2021-10-05', '2021-10-06', '2021-10-07', '2021-10-08', '2021-10-09', '2021-10-10', '2021-10-11', '2021-10-12', '2021-10-13', '2021-10-14', '2021-10-15', '2021-10-16', '2021-10-17', '2021-10-18', '2021-10-19', '2021-10-20', '2021-10-21', '2021-10-22', '2021-10-23', '2021-10-24', '2021-10-25', '2021-10-26', '2021-10-27', '2021-10-28', '2021-10-29', '2021-10-30', '2021-10-31']
NYTNewCases = array([0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,2,0,0,0,3,0,0,1,0,4,1,5,13,12,19,34,59,65,53,60,101,60,102,157,99,83,94,96,134,169,111,75,68,143,53,158,86,65,43,97,133,155,110,110,56,114,123,144,160,133,135,115,127,167,259,204,283,185,178,188,228,91,365,273,63,181,247,129,265,242,283,159,109,160,166,199,132,229,104,111,71,190,280,242,433,319,124,643,563,169,932,719,782,520,408,953,938,957,905,901,748,1562,1563,1713,2351,2345,2019,1616,2478,1184,2276,2614,2867,2784,248,3724,3759,2622,3403,1998,2621,2828,2820,2752,3234,3339,2334,1766,1083,3335,2423,2596,3086,2123,1719,1237,2500,1591,1768,2606,2842,1476,1404,1552,1767,1859,2556,2165,1058,755,832,1286,966,906,650,530,318,551,443,641,633,532,546,244,554,480,492,408,561,140,199,349,159,458,335,429,262,113,342,416,817,481,461,134,56,38,387,303,330,419,228,86,292,315,931,986,-372,289,164,330,120,215,211,244,257,190,417,133,461,292,434,243,229,570,474,585,485,581,450,348,483,570,652,499,655,545,605,698,708,619,649,638,943,731,645,568,815,1013,1196,1073,519,992,620,1382,1411,1591,1215,1089,1561,1005,676,2976,1871,1052,885,3106,2168,3030,2311,2885,936,2364,3390,2297,1497,2666,2525,2000,382,8755,3371,3618,4869,3517,564,394,10321,2800,3606,5233,5031,536,7777,5250,2584,5275,2497,4764,508,5033,6490,3268,5198,929,4161,395,7024,5038,3296,5639,3443,6174,12372,3450,8583,5105,6338,10012,6597,3391,6593,8644,3061,7915,6049,5207,1542,3865,4540,6274,6716,4901,5248,2795,3854,5321,4677,4073,3304,3828,1112,2635,3310,2102,3076,2233,2267,162,1767,4079,1276,2071,856,1252,849,883,765,1261,983,1370,1538,614,1150,1216,787,797,936,643,549,779,840,947,590,1237,1157,1036,605,453,666,1727,1121,53,445,429,383,337,-62,260,279,281,353,418,484,45,388,575,352,385,413,376,226,566,461,313,446,442,571,393,442,603,329,511,490,234,291,733,268,561,516,573,457,349,638,610,520,597,540,468,706,706,870,230,510,541,621,400,690,723,364,517,533,415,403,777,378,426,375,508,466,502,350,574,273,346,424,531,603,710,436,294,264,4,782,266,273,301,568,287,371,448,262,289,472,272,292,346,339,176,254,585,337,301,334,334,473,457,283,323,301,401,578,303,337,516,388,1,719,283,363,683,577,644,85,192,1489,735,943,845,795,754,873,784,810,1172,1188,1205,1149,1155,1047,1281,1559,1502,1990,1501,1618,1624,1833,2102,2142,2084,1814,1995,1600,2196,2336,2607,2320,1953,2087,1736,2608,2344,2338,2627,2036,2114,2440,2439,2852,2651,1671,2455,447,3696,2472,2933,2377,2397,1815,1446,1766,1739,2150,2388,2091,1681,1941,1802,1877,2130,1940,1787,1439,1448,1442,2262,2178,2109,2050,1452,888,2396,1946,2690,2063,1679,1476,1629,1545,1710,1875,1637,1510,1115,1534,1679,1681,1608,1625,1287,1288,1287,1342,1512,2100,1780,1008,73,314,4183,1810,2415,1880,1956])
fig = plt.figure()
fig.set_size_inches(8,5.5)
ax = plt.subplot(1,1,1)
ax2 = ax.twinx()
qtl = load('inferences/Phoenix29-n5-percentiles.npy')
alph = loadmat('prevalence.mat')['alpha_prev_phoenix']
delt = loadmat('prevalence.mat')['delta_prev_phoenix']
tSpan = linspace(0, shape(qtl)[1]-1, shape(qtl)[1])
tSpanExp = linspace(0, len(NYTNewCases)-1, len(NYTNewCases))
ax.scatter(tSpanExp, NYTNewCases, 60, marker='+', linewidth = 1, color = 'black', zorder = 500)
qtlMark = [0.95,0.5,0.2]
ax.fill_between(tSpan, qtl[1,:], qtl[19,:], facecolor = colors1[4], zorder = 0, label=str(int(qtlMark[0]*100))+' % posterior')
ax.fill_between(tSpan, qtl[5,:], qtl[15,:], facecolor = colors1[6], zorder = 1, label=str(int(qtlMark[1]*100))+' % posterior')
ax.fill_between(tSpan, qtl[8,:], qtl[12,:], facecolor = colors1[8], zorder = 2, label=str(int(qtlMark[2]*100))+' % posterior')
ax.set_xticks(linspace(0,len(dateIndex)-1,16))
ax.set_xticklabels(['01-21-20','03-04-20','04-16-20','05-29-20','07-12-20','08-24-20','10-06-20','11-18-20','01-01-21','02-13-21','03-28-21','05-10-21','06-23-21','08-05-21','09-17-21','10-31-21'],rotation=90)
ax.set_xlim((0,len(dateIndex)-1))
ax.set_ylabel('Confirmed case counts')
ax.set_title('Phoenix')
ax.set_ylim((0,13200))
ax2.set_ylim((0,1))
ax2.set_ylabel('Prevalence')
ax2.fill_between(linspace(374,650,40),alph.reshape(40),facecolor = 'yellow', zorder = 0, alpha = 0.5, label='Prevalence of Alpha variant')
ax2.fill_between(linspace(374,650,40),delt.reshape(40),facecolor = 'green', zorder = 0, alpha = 0.5, label='Prevalence of Delta variant')
fig.tight_layout()
fig.legend(loc = (0.24,0.6), frameon=False, fontsize=8)
plt.savefig('Fig5A.pdf')
plt.close()