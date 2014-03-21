source('vash.R')
library(limma)

# Simulate data
# True prior on s_j: Mixture IG
# Unimodal: 'variance': variances~Mix IG with common mode c
#           'precision': precisions~Mix Gamma with common mode c
# n: df of varhat
simvar=function(N,n,pi,alpha,c,unimodal){
  K=length(pi)
  ind=t(rmultinom(N, 1, pi))
  if(unimodal=='variance'){
    modalpha=alpha+1
    truevar=1/rgamma(N*K,shape=rep(alpha,each=N),rate=rep((alpha+1)*c,each=N))
    truevar=apply(ind*matrix(truevar,ncol=K),1,sum)
  }else if(unimodal=='precision'){
    modalpha=alpha-1
    truevar=1/rgamma(N*K,shape=rep(alpha,each=N),rate=rep((alpha-1)*c,each=N))
    truevar=apply(ind*matrix(truevar,ncol=K),1,sum)
  }
  
  Y=matrix(rnorm(N*(n+1),mean=0,sd=rep(sqrt(truevar),n+1)),ncol=n+1)
  varhat=apply(Y,1,var)
  #varhat=rgamma(N,shape=n/2,rate=n/(2*truevar))
  return(list(truevar=truevar,varhat=varhat,Y=Y,ind=ind,n=n,alpha=alpha,modalpha=modalpha,c=c,pi=pi))
}

trueprior=function(dataobj,xmax){
  K=length(dataobj$alpha)
  alpha=dataobj$alpha
  modalpha=dataobj$modalpha
  c=dataobj$c
  pi=dataobj$pi
  xgrid=seq(0.0001,xmax,by=0.01)
  densmat=matrix(dgamma(rep(1/xgrid,K),shape=rep(alpha,each=length(xgrid)),
                        rate=rep(modalpha*c,each=length(xgrid))),ncol=K)
  prior.var=apply(outer(rep(1,length(xgrid)),pi)*densmat,1,sum)*(1/xgrid^2)
  densmat=matrix(dgamma(rep(xgrid,K),shape=rep(alpha,each=length(xgrid)),
                        rate=rep(modalpha*c,each=length(xgrid))),ncol=K)
  prior.prec=apply(outer(rep(1,length(xgrid)),pi)*densmat,1,sum)
  return(list(prior.var=prior.var,prior.prec=prior.prec,xgrid=xgrid))
}

limmastat=function(limmaobj,dataobj,xmax){
  c=limmaobj$s2.prior
  varhat=dataobj$varhat
  alpha.vec=limmaobj$df.prior/2
  modalpha.vec=alpha.vec
  n=dataobj$n
  N=length(varhat)
  pi=1
  K=1
  pimat = exp(outer(rep(1,N),log(pi))+n/2*log(n/2)-lgamma(n/2)
              +(n/2-1)*outer(log(varhat),rep(1,K))
              +outer(rep(1,N),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+n/2))
              -outer(rep(1,N),alpha.vec+n/2)*log(outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  #classprob=pimat/rowSums(pimat)
  logl = sum(log(rowSums(pimat)))
  
  xgrid=seq(0.0001,xmax,by=0.01)
  prior.var=dgamma(xgrid,shape=limmaobj$df.prior/2,rate=limmaobj$df.prior/(2*limmaobj$s2.prior))
  prior.prec=dgamma(1/xgrid,shape=limmaobj$df.prior/2,rate=limmaobj$df.prior/(2*limmaobj$s2.prior))*(1/xgrid^2)
  return(list(loglik=logl,prior.var=prior.var,prior.prec=prior.prec,xgrid=xgrid))  
}
############################
d.1v=simvar(N=10000,n=5,pi=c(0.5,0.5),alpha=c(4,34),c=1,unimodal='variance')
tprior.1v=trueprior(d.1v,3)

vash.1v.v=vash(d.1v$varhat,n=5,unimodal='variance')
vash.1v.vsing=vash(d.1v$varhat,n=5,unimodal='variance',singlecomp=TRUE)
vash.1v.p=vash(d.1v$varhat,n=5,unimodal='precision')

dens.1v.v=vashEBprior(vash.1v.v,3)
plot(dens.1v.v$xgrid,dens.1v.v$EBprior.var,type='l')
#plot(dens.1v.v$xgrid,dens.1v.v$EBprior.prec,type='l')
for (i in 1:dim(dens.1v.v$EBprior.var.sep)[2]){
  lines(dens.1v.v$xgrid,dens.1v.v$EBprior.var.sep[,i],lty=2,col='green')
}

dens.1v.p=vashEBprior(vash.1v.p,3)
#plot(dens.1v.p$xgrid,dens.1v.p$EBprior.var,type='l')
plot(dens.1v.p$xgrid,dens.1v.p$EBprior.prec,type='l')
for (i in 1:dim(dens.1v.p$EBprior.prec.sep)[2]){
  lines(dens.1v.p$xgrid,dens.1v.p$EBprior.prec.sep[,i],lty=2,col='green')
}

# compare with limma
fit.1v <- lmFit(d.1v$Y, rep(1,d.1v$n+1))
fit.1v <- eBayes(fit.1v)
#fit.1v$df.prior/2
#fit.1v$s2.prior
fit.1v.stat=limmastat(fit.1v,d.1v,3)

plot(tprior.1v$xgrid,tprior.1v$prior.var,type='l')
lines(dens.1v.v$xgrid,dens.1v.v$EBprior.var,col='blue')
lines(fit.1v.stat$xgrid,fit.1v.stat$prior.var,col=2)

# True var vs estimated var (vash vs limma)
plot(log2(d.1v$truevar),log2(vash.1v.v$PosteriorMean))
points(log2(d.1v$truevar),log2(fit.1v$s2.post),col=2)
abline(0,1)

mean(abs(log2(d.1v$truevar)-log2(vash.1v.v$PosteriorMean)))
mean(abs(log2(d.1v$truevar)-log2(fit.1v$s2.post)))

# log-likelihood
max(vash.1v.v$pifit$loglik.final,na.rm=TRUE)
fit.1v.stat$loglik

#############
d.1p=simvar(N=10000,n=5,pi=c(0.5,0.5),alpha=c(4,34),c=1,unimodal='precision')
tprior.1p=trueprior(d.1p,3)

vash.1p.v=vash(d.1p$varhat,n=5,unimodal='variance')
vash.1p.p=vash(d.1p$varhat,n=5,unimodal='precision')
dens.1p.v=vashEBprior(vash.1p.v,3)
dens.1p.p=vashEBprior(vash.1p.p,3)

# compare with limma
fit.1p <- lmFit(d.1p$Y, rep(1,d.1p$n+1))
fit.1p <- eBayes(fit.1p)
#fit.1p$df.prior/2
#fit.1p$s2.prior
fit.1p.stat=limmastat(fit.1p,d.1p,3)

plot(tprior.1p$xgrid,tprior.1p$prior.var,type='l')
lines(dens.1p.v$xgrid,dens.1p.v$EBprior.var,col='blue')
lines(fit.1p.stat$xgrid,fit.1p.stat$prior.var,col=2)

plot(log2(d.1p$truevar),log2(vash.1p.p$PosteriorMean),ylim=c(-4,4))
points(log2(d.1p$truevar),log2(fit.1p$s2.post),col=2)
abline(0,1)

mean(abs(log2(d.1p$truevar)-log2(vash.1p.p$PosteriorMean)))
mean(abs(log2(d.1p$truevar)-log2(fit.1p$s2.post)))

max(vash.1p.p$pifit$loglik.final,na.rm=TRUE)
fit.1p.stat$loglik

###############
d.2p=simvar(N=10000,n=5,pi=c(0.5,0.5),alpha=c(18,34),c=3,unimodal='precision')
tprior.2p=trueprior(d.2p,2*d.2p$c)

vash.2p.v=vash(d.2p$varhat,n=5,unimodal='variance')
vash.2p.p=vash(d.2p$varhat,n=5,unimodal='precision')
dens.2p.v=vashEBprior(vash.2p.v,2*d.2p$c)
dens.2p.p=vashEBprior(vash.2p.p,2*d.2p$c)

# compare with limma
fit.2p <- lmFit(d.2p$Y, rep(1,d.2p$n+1))
fit.2p <- eBayes(fit.2p)
#fit.2p$df.prior/2
#fit.2p$s2.prior
fit.2p.stat=limmastat(fit.2p,d.2p,2*d.2p$c)

plot(tprior.2p$xgrid,tprior.2p$prior.var,type='l')
lines(dens.2p.p$xgrid,dens.2p.p$EBprior.var,col='blue')
lines(dens.2p.v$xgrid,dens.2p.v$EBprior.var,col='green')
lines(fit.2p.stat$xgrid,fit.2p.stat$prior.var,col=2)

plot(log2(d.3p$truevar),log2(vash.3p.v$PosteriorMean),ylim=c(-4,4))
points(log2(d.3p$truevar),log2(fit.3p$s2.post),col=2)
abline(0,1)

mean(abs(log2(d.2p$truevar)-log2(vash.2p.p$PosteriorMean)))
mean(abs(log2(d.2p$truevar)-log2(fit.2p$s2.post)))

max(vash.2p.p$pifit$loglik.final,na.rm=TRUE)
fit.2p.stat$loglik

###########################
d.3p=simvar(N=10000,n=5,pi=c(0.2,0.8),alpha=c(4,34),c=3,unimodal='precision')
tprior.3p=trueprior(d.3p,2*d.3p$c)

vash.3p.v=vash(d.3p$varhat,n=5,unimodal='variance')
vash.3p.p=vash(d.3p$varhat,n=5,unimodal='precision')
dens.3p.v=vashEBprior(vash.3p.v,2*d.3p$c)
dens.3p.p=vashEBprior(vash.3p.p,2*d.3p$c)

# compare with limma
fit.3p <- lmFit(d.3p$Y, rep(1,d.3p$n+1))
fit.3p <- eBayes(fit.3p)
#fit.3p$df.prior/2
#fit.3p$s2.prior
fit.3p.stat=limmastat(fit.3p,d.3p,2*d.3p$c)

plot(tprior.3p$xgrid,tprior.3p$prior.var,type='l')
lines(dens.3p.v$xgrid,dens.3p.v$EBprior.var,col='blue')
lines(fit.3p.stat$xgrid,fit.3p.stat$prior.var,col=2)

plot(log2(d.3p$truevar),log2(vash.3p.p$PosteriorMean))
points(log2(d.3p$truevar),log2(fit.3p$s2.post),col=2)
abline(0,1)

mean(abs(log2(d.3p$truevar)-log2(vash.3p.p$PosteriorMean)))
mean(abs(log2(d.3p$truevar)-log2(fit.3p$s2.post)))

max(vash.3p.p$pifit$loglik.final,na.rm=TRUE)
fit.3p.stat$loglik

###########################
d.4p=simvar(N=10000,n=5,pi=c(0.8,0.2),alpha=c(4,34),c=3,unimodal='precision')
tprior.4p=trueprior(d.4p,2*d.4p$c)

vash.4p.v=vash(d.4p$varhat,n=5,unimodal='variance')
vash.4p.p=vash(d.4p$varhat,n=5,unimodal='precision')
dens.4p.v=vashEBprior(vash.4p.v,2*d.4p$c)
dens.4p.p=vashEBprior(vash.4p.p,2*d.4p$c)

# compare with limma
fit.4p <- lmFit(d.4p$Y, rep(1,d.4p$n+1))
fit.4p <- eBayes(fit.4p)
#fit.4p$df.prior/2
#fit.4p$s2.prior
fit.4p.stat=limmastat(fit.4p,d.4p,2*d.4p$c)

plot(tprior.4p$xgrid,tprior.4p$prior.var,type='l')
lines(dens.4p.v$xgrid,dens.4p.v$EBprior.var,col='blue')
lines(fit.4p.stat$xgrid,fit.4p.stat$prior.var,col=2)

plot(log2(d.4p$truevar),log2(vash.4p.v$PosteriorMean),ylim=c(-4,4))
points(log2(d.4p$truevar),log2(fit.4p$s2.post),col=2)
abline(0,1)

mean(abs(log2(d.4p$truevar)-log2(vash.4p.p$PosteriorMean)))
mean(abs(log2(d.4p$truevar)-log2(fit.4p$s2.post)))

max(vash.4p.p$pifit$loglik.final,na.rm=TRUE)
fit.4p.stat$loglik

#####################
d.5p=simvar(N=10000,n=50,pi=c(0.2,0.8),alpha=c(4,34),c=3,unimodal='precision')
tprior.5p=trueprior(d.5p,2*d.5p$c)

#vash.5p.v=vash(d.5p$varhat,n=50,unimodal='variance')
vash.5p.p=vash(d.5p$varhat,n=50,unimodal='precision')
#dens.5p.v=vashEBprior(vash.5p.v,2*d.5p$c)
dens.5p.p=vashEBprior(vash.5p.p,2*d.5p$c)

# compare with limma
fit.5p <- lmFit(d.5p$Y, rep(1,d.5p$n+1))
fit.5p <- eBayes(fit.5p)
#fit.5p$df.prior/2
#fit.5p$s2.prior
fit.5p.stat=limmastat(fit.5p,d.5p,2*d.5p$c)

plot(tprior.5p$xgrid,tprior.5p$prior.var,type='l')
lines(dens.5p.v$xgrid,dens.5p.v$EBprior.var,col='blue')
lines(fit.5p.stat$xgrid,fit.5p.stat$prior.var,col=2)

plot(log2(d.5p$truevar),log2(vash.5p.v$PosteriorMean))
points(log2(d.5p$truevar),log2(fit.5p$s2.post),col=2)
abline(0,1)

mean(abs(log2(d.5p$truevar)-log2(vash.5p.p$PosteriorMean)))
mean(abs(log2(d.5p$truevar)-log2(fit.5p$s2.post)))

max(vash.5p.p$pifit$loglik.final,na.rm=TRUE)
fit.5p.stat$loglik

#####################
d.6p=simvar(N=10000,n=5,pi=rep(0.1,10),alpha=seq(3,48,by=5),c=2,unimodal='precision')
tprior.6p=trueprior(d.6p,2*d.6p$c)

#vash.6p.v=vash(d.6p$varhat,n=5,unimodal='variance')
vash.6p.p=vash(d.6p$varhat,n=5,unimodal='precision')
#dens.6p.v=vashEBprior(vash.6p.v,2*d.6p$c)
dens.6p.p=vashEBprior(vash.6p.p,2*d.6p$c)

# compare with limma
fit.6p <- lmFit(d.6p$Y, rep(1,d.6p$n+1))
fit.6p <- eBayes(fit.6p)
#fit.6p$df.prior/2
#fit.6p$s2.prior
fit.6p.stat=limmastat(fit.6p,d.6p,2*d.6p$c)

plot(tprior.6p$xgrid,tprior.6p$prior.var,type='l')
lines(dens.6p.p$xgrid,dens.6p.p$EBprior.var,col='blue')
lines(fit.6p.stat$xgrid,fit.6p.stat$prior.var,col=2)

plot(log2(d.6p$truevar),log2(vash.6p.p$PosteriorMean),ylim=c(-4,4))
points(log2(d.6p$truevar),log2(fit.6p$s2.post),col=2)
abline(0,1)

mean(abs(log2(d.6p$truevar)-log2(vash.6p.p$PosteriorMean)))
mean(abs(log2(d.6p$truevar)-log2(fit.6p$s2.post)))

max(vash.6p.p$pifit$loglik.final,na.rm=TRUE)
fit.6p.stat$loglik

#####################
d.7p=simvar(N=10000,n=20,pi=c(rep(0.1,9),rep(0.01,10)),alpha=seq(4,40,by=2),c=3,unimodal='precision')
tprior.7p=trueprior(d.7p,2*d.7p$c)

#vash.7p.v=vash(d.7p$varhat,n=20,unimodal='variance')
vash.7p.p=vash(d.7p$varhat,n=20,unimodal='precision')
#dens.7p.v=vashEBprior(vash.7p.v,2*d.7p$c)
dens.7p.p=vashEBprior(vash.7p.p,2*d.7p$c)

# compare with limma
fit.7p <- lmFit(d.7p$Y, rep(1,d.7p$n+1))
fit.7p <- eBayes(fit.7p)
#fit.7p$df.prior/2
#fit.7p$s2.prior
fit.7p.stat=limmastat(fit.7p,d.7p,2*d.7p$c)

plot(tprior.7p$xgrid,tprior.7p$prior.var,type='l')
lines(dens.7p.p$xgrid,dens.7p.p$EBprior.var,col='blue')
lines(fit.7p.stat$xgrid,fit.7p.stat$prior.var,col=2)

plot(log2(d.7p$truevar),log2(vash.7p.p$PosteriorMean),ylim=c(-4,4))
points(log2(d.7p$truevar),log2(fit.7p$s2.post),col=2)
abline(0,1)

mean(abs(log2(d.7p$truevar)-log2(vash.7p.p$PosteriorMean)))
mean(abs(log2(d.7p$truevar)-log2(fit.7p$s2.post)))

max(vash.7p.p$pifit$loglik.final,na.rm=TRUE)
fit.7p.stat$loglik
