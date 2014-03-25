install.packages("../../rcpp-package/ashr_0.1.tar.gz",repos=NULL,type="source")
library(ashr)

# simulate n beta-hat values, nnull under the null with altmean and altsd
# being the mean and sd of beta under the alternative
simdata = function(n, nnull, altmean, altsd, betasd) {
        null = c(rep(1, nnull), rep(0, n - nnull))
            beta = c(rep(0, nnull), rnorm(n - nnull, altmean, altsd))
            betahat = rnorm(n, beta, betasd)
            return(list(null = null, beta = beta, betahat = betahat, betasd = betasd))
    }

set.seed(100)
ss = simdata(10000, 8000, 0, 2, 1)

# do profiling
cxx=TRUE
outfile="ash.cxx.out"
Rprof(outfile)
system.time(beta.ash.cxx <- ash(ss$betahat, ss$betasd, prior = "uniform", cxx =cxx))
Rprof(NULL)

cxx=FALSE
outfile="ash.out"
Rprof(outfile)
system.time(beta.ash <- ash(ss$betahat, ss$betasd, prior = "uniform", cxx =cxx))
Rprof(NULL)


# check that return identical values (up to machine error)
all.equal(beta.ash$PosteriorMean, beta.ash.cxx$PosteriorMean)
all.equal(beta.ash$qvalue, beta.ash.cxx$qvalue)
all.equal(beta.ash$PositiveProb, beta.ash.cxx$PositiveProb)

#to see profiling output run this in the terminal
#R CMD Rprof ash.out > ash.out.txt
#R CMD Rprof ash.cxx.out > ash.cxx.out.txt

# repeat using halfuniform
cxx=TRUE
outfile="ash.cxx.out"
Rprof(outfile)
system.time(beta.ash.cxx <- ash(ss$betahat, ss$betasd, prior = "uniform", mixcompdist="halfuniform",cxx =cxx))
Rprof(NULL)

cxx=FALSE
outfile="ash.out"
Rprof(outfile)
system.time(beta.ash <- ash(ss$betahat, ss$betasd, prior = "uniform", mixcompdist="halfuniform", cxx =cxx))
Rprof(NULL)

# check that return identical values (up to machine error)
all.equal(beta.ash$PosteriorMean, beta.ash.cxx$PosteriorMean)
all.equal(beta.ash$qvalue, beta.ash.cxx$qvalue)
all.equal(beta.ash$PositiveProb, beta.ash.cxx$PositiveProb)

# check that convergence warning is produced
ash(ss$betahat, ss$betasd, prior = "uniform", mixcompdist="halfuniform",cxx =TRUE,maxiter=2)
ash(ss$betahat, ss$betasd, prior = "uniform", mixcompdist="halfuniform",cxx =FALSE,maxiter=2)
ash(ss$betahat, ss$betasd, prior = "uniform", mixcompdist="halfuniform",VB=TRUE,cxx =FALSE,maxiter=2)

