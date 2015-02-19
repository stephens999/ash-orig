#R code to run some simulations

require(ashr)
require(qvalue)
require(fdrtool)
require(mixfdr)
require(locfdr)
require(ggplot2)
#bsd gives standard deviation of beta
#pi0 is simulated to be uniform on [minpi0,1]
basicsim=function(mixsd,mixpi_alt,bsd=1,minpi0=0,seedval = 100,nsamp=1000,niter=50){  
  set.seed(seedval)  
  beta =list()
  betahatsd=list()
  betahat = list()
  zscore = list()
  pval = list()
  betahat.ash.n = list()
  betahat.ash.u = list()
  betahat.ash.ur = list()
  betahat.ash.npm = list()
  betahat.ash.true = list()
  betahat.qval = list()
  betahat.fdrtool = list()
  betahat.locfdr = list()
  betahat.mixfdr = list()
  pi0 = rep(0,niter)
  for(i in 1:niter){
    pi0[i]=runif(1,minpi0,1)
    mixpi = c(pi0[i],(1-pi0[i])*mixpi_alt)
    sd = sample(mixsd,nsamp,prob=mixpi,replace=TRUE )
    beta[[i]] = rnorm(nsamp,0,sd)
    betahatsd[[i]] = bsd
    betahat[[i]] = beta[[i]]+rnorm(nsamp,0,betahatsd[[i]])
    zscore[[i]] = betahat[[i]]/betahatsd[[i]]
    pval[[i]] = pchisq(zscore[[i]]^2,df=1,lower.tail=F)
    betahat.ash.n[[i]] = ash(betahat[[i]],betahatsd[[i]],pointmass=TRUE,prior="nullbiased",gridmult=2)
    betahat.ash.u[[i]] = ash(betahat[[i]],betahatsd[[i]],pointmass=TRUE,prior="uniform",gridmult=2)
    betahat.ash.ur[[i]] = ash(betahat[[i]],betahatsd[[i]],pointmass=TRUE,prior="uniform",randomstart=TRUE,gridmult=2)
    betahat.ash.npm[[i]] = ash(betahat[[i]],betahatsd[[i]],pointmass=FALSE,prior="uniform",gridmult=2)
    betahat.ash.true[[i]] = ash(betahat[[i]],betahatsd[[i]],g=normalmix(mixpi,rep(0,length(mixpi)),mixsd),control=list(maxiter=0))
    
    cat("applying q value\n")
    betahat.qval[[i]] = qvalue(pval[[i]])
    cat("applying fdrtool\n")
    betahat.fdrtool[[i]] = fdrtool(pval[[i]],statistic="pvalue",plot=FALSE)
    cat("applying locfdr\n")
    betahat.locfdr[[i]] = try(locfdr(zscore[[i]],nulltype=0,plot=0))
    cat("applying mixfdr\n")
    betahat.mixfdr[[i]] = mixFdr(zscore[[i]],noiseSD=1,theonull=TRUE,plot=FALSE)
  }
  return(list(beta =beta,
              betahatsd=betahatsd,
              betahat = betahat,
              zscore = zscore,
              pval = pval,
              betahat.ash.n = betahat.ash.n,
              betahat.ash.u = betahat.ash.u,
              betahat.ash.npm = betahat.ash.npm,
              betahat.ash.true=betahat.ash.true,
              betahat.qval = betahat.qval,
              betahat.fdrtool = betahat.fdrtool,
              betahat.locfdr = betahat.locfdr,
              betahat.mixfdr = betahat.mixfdr,
              pi0=pi0))
}    

mixsd = c(0,0.25,0.5,1,2)
mixpi_alt = c(0.4,0.2,0.2,0.2) #mixture proportions under the alternative

set.seed(100)
#these are the simulations for scenarios 1a, 1b and 2 from the paper
simres1a = basicsim(mixsd,mixpi_alt,niter=200,nsamp=1000)
simres1b = basicsim(mixsd,mixpi_alt,niter=200,nsamp=10000)
simres2= basicsim(c(0,4),c(1))

#these do a situation where precision varies across measurements
simres3=basicsim(c(0,2),c(1),bsd=c(rep(1,500),rep(10,500)))
simres4=basicsim(c(0,2),c(1),bsd=c(rep(1,500),rep(10,500)),minpi0=0.9,seed=200)

save.image(file="sim1.RData")
