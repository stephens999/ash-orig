

load("debug_logLR_inf.Robj")

# load library ashr

A1 = ash(zdat.rate[3],zdat.rate[4], prior=prior, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR = TRUE)
A2 = fast.ash(zdat.rate[3],zdat.rate[4], prior="uniform", pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB, onlylogLR = TRUE)


A1$logLR
#[1] Inf
A2$logLR
#[1] Inf

