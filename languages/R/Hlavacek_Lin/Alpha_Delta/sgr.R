setwd("~/Documents/bmab")
library(reticulate)
library(stableGR)
np <- import("numpy")
x1 <- np$load('inferences/Dallas4-n5-parLog.npy')[220000:340000,]
x2 <- np$load('inferences/Houston2-n4-parLog.npy')[200000:380000,]
x3 <- np$load('inferences/NYC0-n4-parLog.npy')[1:40000,]
x4 <- np$load('inferences/Phoenix29-n5-parLog.npy')[1:100000,]
print(stable.GR(x1)$mpsrf)
print(stable.GR(x1)$n.eff)
print(stable.GR(x2)$mpsrf)
print(stable.GR(x2)$n.eff)
print(stable.GR(x3)$mpsrf)
print(stable.GR(x3)$n.eff)
print(stable.GR(x4)$mpsrf)
print(stable.GR(x4)$n.eff)