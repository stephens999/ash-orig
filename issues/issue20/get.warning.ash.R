



library("ashr")

load("get.warning.ash.Robj")
 
zdat.ash = ash(zdat[3,ind], zdat[4,ind], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR=onlylogLR)
#Warning message:
#In EMest(betahat[completeobs], lambda1 * sebetahat[completeobs] +  :
#  EM algorithm in function mixEM failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.

zdat.ash$logLR
#[1] 15475.21
zdat.ash$fitted.g$pi
# [1] 5.042176e-03 5.830010e-08 1.946892e-07 1.395479e-05 6.376346e-01
# [6] 1.936746e-02 3.045065e-01 3.343449e-02 5.742198e-07 0.000000e+00
zdat.ash$fitted.g$sd
# [1] 0.00000000 0.01220149 0.02440298 0.04880597 0.09761193 0.19522387
# [7] 0.39044773 0.78089547 1.56179094 3.12358188
