#simulate from mixture of normals
#compare normal, uniform and half-uniform mixtures
#bsd gives standard deviation of beta
require(ashr)
require(mixfdr)

ashsim=function(mixmean,mixsd,mixpi,bsd=1,seedval = 100,nsamp=1000,niter=50){  
  set.seed(seedval)  
  beta =list()
  betahatsd=list()
  betahat = list()
  fit.ash.n = list()
  fit.ash.u = list()
  fit.ash.hu = list()
  fit.ash.hu.vb=list()
  fit.ash.true = list()
  fit.mixfdr= list()
  fit.ash.fdr.n = list()
  fit.ash.fdr.u = list()
  fit.ash.fdr.hu = list()
  
  
  fit.mixfdr.enull= list()
  fit.mixfdr.J10= list()
  fit.mixfdr.J10P0= list()
  fit.mixfdr.J100= list()
  k = length(mixmean)
  for(i in 1:niter){
    comp = sample(1:k,nsamp,prob=mixpi,replace=TRUE)
    sd = mixsd[comp]
    mean = mixmean[comp]
    beta[[i]] = rnorm(nsamp,mean,sd)
    betahatsd[[i]] = bsd
    betahat[[i]] = beta[[i]]+rnorm(nsamp,0,betahatsd[[i]])
    fit.ash.n[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="normal",method="shrink")
    fit.ash.u[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="uniform", method="shrink")
    fit.ash.hu[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="halfuniform", method="shrink")
    fit.ash.hu.vb[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="halfuniform",method="shrink",VB=TRUE)
    fit.ash.true[[i]] = ash(betahat[[i]],betahatsd[[i]],g=normalmix(mixpi,mixmean,mixsd),control=list(maxiter=0))
    fit.mixfdr[[i]] = mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=TRUE,plot=FALSE)
    
    fit.ash.fdr.n[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="normal",method="fdr")
    fit.ash.fdr.u[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="uniform", method="fdr")
    fit.ash.fdr.hu[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="halfuniform", method="fdr")
    fit.mixfdr.enull[[i]] = mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=FALSE,plot=FALSE)
    fit.mixfdr.J10[[i]] = mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=TRUE,plot=FALSE,J=10)
    #fit.mixfdr.J100[[i]] = mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=TRUE,plot=FALSE,J=100)
    fit.mixfdr.J10P0[[i]] = try(mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=TRUE,plot=FALSE,J=10,P=0))
  }
  return(list(beta =beta,
              betahatsd=betahatsd,
              betahat = betahat,
              fit.ash.n = fit.ash.n,
              fit.ash.u = fit.ash.u,
              fit.ash.hu = fit.ash.hu,
              fit.ash.hu.vb=fit.ash.hu.vb,
              fit.mixfdr = fit.mixfdr,
              fit.ash.fdr.n = fit.ash.fdr.n,
              fit.ash.true=fit.ash.true,
              # fit.ash.fdr.u = fit.ash.fdr.u,
              # fit.ash.fdr.hu = fit.ash.fdr.hu,
              fit.mixfdr.enull = fit.mixfdr.enull,
              fit.mixfdr.J10 = fit.mixfdr.J10,
              #fit.mixfdr.J100= fit.mixfdr.J100,
              fit.mixfdr.J10P0=fit.mixfdr.J10P0))
} 

set.seed(111)
simA= ashsim(c(0,0,0),c(1,1,2),c(1/3,1/3,1/3),niter=100,nsamp=1000)
simB= ashsim(c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7),rep(1/7,7),niter=100,nsamp=1000)
simC= ashsim(c(-2,-1,0,1),c(2,1.5,1,1),c(1/4,1/4,1/3,1/6),niter=100,nsamp=1000)
save.image(file="simABC.RData")
