

load("../issues/issue10/debug_logLR_inf.Robj")

# load library ashr

A1 = ash(zdat.rate[3],zdat.rate[4])
A2 = fast.ash(zdat.rate[3],zdat.rate[4],onlylogLR=TRUE)

A1$logLR
#[1] Inf
A2$logLR
#[1] Inf

set.seed(100)
A1 = ash(c(zdat.rate[3],rnorm(1000)),c(zdat.rate[4],rep(1,1000)))
A1$logLR
