library(Sim.DiffProc)
library(plyr)
library(R.matlab)

Mode <- function(x) 
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

set.seed(1234, kind = "L'Ecuyer-CMRG")
St <- expression(0)
k <- 8
arr_b1 <- seq(0.01,0.99,length.out=k)
arr_b2 <- seq(0.01,0.99,length.out=k)
arr_b3 <- seq(0.01,0.99,length.out=k)
arr_d <- seq(0.01,0.99,length.out=k)
arr_sa <- c(0.1,0.5,1)
pars <- expand.grid(b1=arr_b1,b2=arr_b2,b3=arr_b3,d=arr_d,sa=arr_sa)
adat2 <- matrix(NA,nrow(pars),6)

for (i in 1:nrow(pars))
{
  fx <- expression(pars$b1[i]*x+2*x^2-x-x^3+pars$d[i]*(y+z-2*x),pars$b2[i]*y+2*y^2-y-y^3+pars$d[i]*(x+z-2*y),pars$b3[i]*z+2*z^2-z-z^3+pars$d[i]*(x+y-2*z))
  gx <- expression(pars$sa[i],pars$sa[i],pars$sa[i])
  mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(x=1+sqrt(pars$b1[i]),y=1+sqrt(pars$b2[i]),z=1+sqrt(pars$b3[i])),M=1000,T=1000)
  fpt3d <- fptsde3d(mod3d, boundary = St)
  outx <- fpt3d$fpt$x
  outy <- fpt3d$fpt$y
  outz <- fpt3d$fpt$z  
  mfptx <- mean(outx)
  mfpty <- mean(outy)
  mfptz <- mean(outz)
  lfptx <- Mode(outx)
  lfpty <- Mode(outy)
  lfptz <- Mode(outz)
  adat2[i,] <- c(mfptx, mfpty, mfptz, lfptx, lfpty, lfptz)
  print(i)
}

writeMat('additive_T2.mat',adat2=adat2)