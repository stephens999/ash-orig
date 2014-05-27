# Y: N by n observation matrix (N tests, n replicate obs for each test)
# Model: Y_{jk}|beta_j,tau_j ~ N(beta_j,1/tau_j)
# Prior: beta_j|tau_j~N(mu_j,1/(c_l*tau_j)), tau_j~lambda_k*Gamma(a_m,a_m), w.p. pi_{klm}
# Mixture constraint: sum(pi_{klm})=1
# Or: sum(pi_k)=1, sum(pi_l)=1, sum(pi_m)=1 for each param
# a_m, lambda_k, c_l, mu_j are known
######################################################
library(qvalue)
qval.from.localfdr = function(localfdr){
  o = order(localfdr)
  qvalue=rep(NA,length(localfdr))
  qvalue[o] = (cumsum(sort(localfdr))/(1:sum(!is.na(localfdr))))
  return(qvalue)
}


# post_pi computes the posterior probability 
# P(beta_j,tau_j belongs to the i'th component | Y)
post_pi = function(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE,pi){
  post.gammab=0.5*(outer(SSE,rep(1,M*L*K))+outer((MEAN-mu)^2,n*d.vec)+outer(rep(1,N),2*b.vec))
  post.pi.mat = exp(outer(rep(1,N),log(pi)+a.vec*log(b.vec)-lgamma(a.vec)+0.5*log(d.vec)+lgamma(a.vec+n/2))-
                      log(post.gammab)*
                      outer(rep(1,N),a.vec+n/2))
  return(list(gammab=post.gammab,pimat=post.pi.mat))
  
#   post.pi.mat = exp(outer(rep(1,N),log(pi)+a.vec*log(b.vec)-lgamma(a.vec)+log(sqrt(d.vec)+lgamma(a.vec+n/2))-
#     log(0.5*(outer(SSE,rep(1,M*L*K))+outer((MEAN-mu)^2,n*d.vec)+outer(rep(1,N),2*b.vec)))*
#     outer(rep(1,N),a.vec+n/2))
 
#   post.pi.mat = outer(rep(1,N),pi*(b.vec)^a.vec/gamma(a.vec)*sqrt(d.vec*(gamma(a.vec+n/2)))/
#     (0.5*(outer(SSE,rep(1,M*L*K))+outer((MEAN-mu)^2,n*d.vec)+outer(rep(1,N),2*b.vec)))^
#     outer(rep(1,N),a.vec+n/2)
}



# indep_post_pi computes the posterior probability 
# P(beta_j,tau_j belongs to the i'th component of precShape/precMulti.compprecPrior| Y)
# groupind: i indicates that the current component correspond to the i'th comp of precShape/precMulti/compprecPrior
indep_post_pi = function(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE, pi.a, pi.lambda, pi.c, prior, groupind){
  pi = rep(pi.a,L*K)*rep(rep(pi.lambda,each=M),L)*rep(pi.c,each=M*K)
  mm=post_pi(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE,pi)
  pi.mat = mm$pimat
  postb=mm$gammab
  norm.pi.mat = normalized_pi(pi.mat)
  indep.post.pi = rep(0,max(groupind))
  temp = colSums(norm.pi.mat)+prior-5/(M*K)
  #temp=colSums(norm.pi.mat)
  for (i in 1:max(groupind)){
    indep.post.pi[i] = sum(temp[groupind==i])
  } 
  indep.post.pi = ifelse(indep.post.pi<0,0,indep.post.pi)
  indep.post.pi=indep.post.pi/sum(indep.post.pi)
  return(list(indep.pi=indep.post.pi,pi.mat=pi.mat,normpi.mat=norm.pi.mat,postb=postb))
}

# Normalize pi to make sum(pi)=1
normalized_pi = function(pi.mat){
  n = dim(pi.mat)[2]
  pi.normalized = pi.mat/outer(rowSums(pi.mat),rep(1,n))  
  return(pi.normalized)
}

# Compute posterior distribution P(tau|Y), P(beta|Y)
# pi: (M*L*K) vector
# mu: N vector
post_distn = function(N,n,M,K,L,a.vec,b.vec,c.vec,mu,MEAN,SSE,pi){ 
  d.vec=1/(n/c.vec+1)
  post.pi = normalized_pi(post_pi(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE,pi)$pimat)

  post.gamma.parama = outer(rep(1,N),a.vec+n/2)
  post.gamma.paramb = 0.5*(outer(SSE,rep(1,M*L*K))+outer((MEAN-mu)^2,n*d.vec)+outer(rep(1,N),2*b.vec))
  
  # Posterior mean: E(tau|Y)
  post.tau = apply(post.gamma.parama/post.gamma.paramb*post.pi,1,sum)/apply(post.pi,1,sum)
  
  # Posterior mean: E(beta|Y)
  post.norm.mean = outer(MEAN,1-d.vec)+outer(mu,d.vec)
  post.norm.prec = n+c.vec # if tau given, outer(tau,n+c.vec)
  post.norm.c = outer(rep(1,N),n+c.vec)
  post.gamma.dens = dgamma(outer(post.tau,rep(1,M*L*K)),shape=post.gamma.parama,rate=post.gamma.paramb)
  post.beta = apply(post.pi*post.norm.mean,1,sum)
  
  return(list(pi=post.pi,tau=post.tau,beta=post.beta,gammaa=post.gamma.parama,gammab=post.gamma.paramb,
              gammadens=post.gamma.dens,normmean=post.norm.mean,normc=post.norm.c,normprec=post.norm.prec))
}

#loglike=function(params){
loglike=function(params,N,n,M,K,L,mu,MEAN,SSE,pi,groupind){
  a.vec = rep(params[1:M],L*K)
  b.vec = rep(rep(params[(M+1):(M+K)],each=M),L)
  #d.vec = rep(params[(M+K+1):(M+K+L)],each=M*K)
  d.vec = rep(c(params[(M+K+1):(M+K+L-1)],1),each=M*K)
  post.gammab=0.5*(outer(SSE,rep(1,M*L*K))+outer((MEAN-mu)^2,n*d.vec)+outer(rep(1,N),2*b.vec))
  pimat=exp(outer(rep(1,N),log(pi)+a.vec*log(b.vec)-lgamma(a.vec)+0.5*log(d.vec)+lgamma(a.vec+n/2))-
    log(post.gammab)*
    outer(rep(1,N),a.vec+n/2))
  #classprob=pimat/rowSums(pimat)
  logl=sum(log(rowSums(pimat)))
  
#   grad=rep(0,K+M+L)
#   for (i in 1:M){
#     idx=(groupind[,1]==i)
#     a=params[i]
#     grad[i]=sum(classprob[,idx]*(outer(rep(1,N),log(b.vec[idx]))-log(post.gammab[,idx])+outer(rep(1,N),rep(digamma(a+n/2)-digamma(a),length(idx)))))
#   }
#   for (i in 1:K){
#     idx=(groupind[,2]==i)
#     b=params[i+M]
#     grad[i+M]=sum(classprob[,idx]*(outer(rep(1,N),a.vec[idx]/b)-outer(rep(1,N),a.vec[idx]+n/2)/post.gammab[,idx]))
#   }
#   for (i in 1:L){
#     idx=(groupind[,3]==i)
#     d=params[i+M+K]
#     grad[i+M+K]=sum(classprob[,idx]*(1/(2*d)-outer(n*(MEAN-mu)^2/2,a.vec[idx]+n/2)/post.gammab[,idx]))    
#   }
#   attr(logl, "gradient") = grad
  return(-logl)
}

gradloglike=function(params,N,n,M,K,L,mu,MEAN,SSE,pi,groupind){
  a.vec = rep(params[1:M],L*K)
  b.vec = rep(rep(params[(M+1):(M+K)],each=M),L)
  d.vec = rep(c(params[(M+K+1):(M+K+L-1)],1),each=M*K)
  #d.vec = rep(params[(M+K+1):(M+K+L)],each=M*K)
  post.gammab=0.5*(outer(SSE,rep(1,M*L*K))+outer((MEAN-mu)^2,n*d.vec)+outer(rep(1,N),2*b.vec))
  pimat=exp(outer(rep(1,N),log(pi)+a.vec*log(b.vec)-lgamma(a.vec)+0.5*log(d.vec)+lgamma(a.vec+n/2))-
              log(post.gammab)*
              outer(rep(1,N),a.vec+n/2))
  classprob=pimat/rowSums(pimat)
  grad=rep(0,K+M+L-1)
  for (i in 1:M){
    idx=(groupind[,1]==i)
    a=params[i]
    grad[i]=-sum(classprob[,idx]*(outer(rep(1,N),log(b.vec[idx]))-log(post.gammab[,idx])+outer(rep(1,N),rep(digamma(a+n/2)-digamma(a),length(idx)))))
  }
  for (i in 1:K){
    idx=(groupind[,2]==i)
    b=params[i+M]
    grad[i+M]=-sum(classprob[,idx]*(outer(rep(1,N),a.vec[idx]/b)-outer(rep(1,N),a.vec[idx]+n/2)/post.gammab[,idx]))
  }
  for (i in 1:(L-1)){
    idx=(groupind[,3]==i)
    d=params[i+M+K]
    grad[i+M+K]=-sum(classprob[,idx]*(1/(2*d)-outer(n*(MEAN-mu)^2/2,a.vec[idx]+n/2)/post.gammab[,idx]))    
  }
  return(grad)
}

# Use EM algorithm to estimate a, lambda
# lik: likelihood matrix
# SGD: use stochastic gradient descent method
a_lambda_c_est=function(classprob,postb,N,n,M,K,L,a.vec,b.vec,d.vec,groupind,SGD,learnrate){
  if(SGD==TRUE){    
    dl_da=sum(classprob*(outer(rep(1,N),log(b.vec))-log(postb)+outer(rep(1,N),digamma(a.vec+n/2)-digamma(a.vec))))
    dl_db=sum(classprob*(outer(rep(1,N),a.vec/b.vec)-outer(rep(1,N),a.vec+n/2)/postb))
    a.vec=pmin(a.vec-learnrate*dl_da,1e-10)
    b.vec=pmin(b.vec-learnrate*dl_db,1e-10) 
    print(dl_da)
    print(dl_db)
  }
  else{
    for(i in 1:max(group.ind[,1])){
      idx=(group.ind[,1]==i)
      fa=function(a) {        
        res=sum(classprob[,idx]*(outer(rep(1,N),log(b.vec[idx]))-log(postb[,idx])+outer(rep(1,N),rep(digamma(a+n/2)-digamma(a),length(idx)))))
        return(res)
      }
      a.vec[idx]=uniroot(fa,c(1e-10,100000))$root   
    }
    for(i in 1:max(group.ind[,2])){
      idx=(group.ind[,2]==i)
      fb=function(b) {
        res=sum(classprob[,idx]*(outer(rep(1,N),a.vec[idx]/b)-outer(rep(1,N),a.vec[idx]+n/2)/postb[,idx]))
        return(res)
      }
      b.vec[idx]=uniroot(fb,c(1e-10,100000))$root   
    }
#     for(i in 1:max(group.ind[,3])){
#       idx=(group.ind[,3]==i)
#       fd=function(d) {
#         res=sum(classprob[,idx]*(1/(2*d)-outer(n*(MEAN-mu)^2/2,a.vec[idx]+n/2)/postb[,idx]))
#         return(res)
#       }
#       optd=uniroot(fd,c(1e-10,1e10))$root
#       if (optd<1){
#         d.vec[idx]=optd
#       }else{
#         d.vec[idx]=1
#       }    
    for(i in 1:(max(group.ind[,3]))){
      idx=(group.ind[,3]==i)
      fc=function(c) {
        res=sum(classprob[,idx]*(n/(2*c*(n+c))-outer(n/(n+c)^2*0.5*n*(MEAN-mu)^2,a.vec[idx]+n/2)/postb[,idx]))
        return(res)
      }
      c.vec[idx]=uniroot(fd,c(1e-10,1e10))$root
    }
  }
  return(list(a.vec=a.vec,b.vec=b.vec,d.vec=d.vec,c.vec=c.vec))
}



# Use EM algorithm to estimate pi from data
# ltol: likelihood convergence tolerance
# maxiter: max number of iterations
# indepprior: whether assume that there are independent priors for a_m, lambda_k and c_l
# i.e. pi_{klm}=pi.lambda_k*pi.c_l*pi.a_m, where sum(pi.a_m)=sum(pi.lambda_k)=sum(pi.c_l)
EMest_pi = function(params,N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE, pi, prior, a.lambda.c.est, SGD,indepprior=TRUE,ltol=0.0001, maxiter=5000, usePointMass){
  group.a = rep(1:M,L*K)
  group.lambda = rep(rep(1:K,each=M),L)
  group.c = rep(1:L,each=M*K)
  groupind=cbind(group.a,group.lambda,group.c)
  
  if(is.null(prior)){ # NEED TESTING
    prior=rep(5/(L*M*K),L*M*K)
    prior[d.vec==max(d.vec)] = 5/(M*K)     
  }else if(prior=="uniform"){
    prior=rep(5/(M*K),M*L*K)
  }
  if(indepprior==FALSE){
    pi.a=NULL; pi.lambda=NULL; pi.c=NULL
    if(is.null(pi)){
      pi=rep(1,K*L*M)/(K*L*M)
      pi[d.vec==max(d.vec)]=L
      pi=pi/sum(pi)
    }
    loglik = rep(NA,maxiter)
    
    mm = post_pi(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE,pi)
    m = mm$pimat
    m.rowsum = rowSums(m)
    loglik[1] = sum(log(m.rowsum))
    classprob = m/m.rowsum
    
    for(i in 2:maxiter){
      if(a.lambda.c.est==TRUE){
        est=nlminb(params,loglike,gradloglike,lower=rep(1e-10,M+K+L),upper=c(rep(100000,M+K),rep(0.95,L)),N=N,n=n,M=M,K=K,L=L,mu=mu,MEAN=MEAN,SSE=SSE,pi=pi,groupind=groupind)
        #est=nlminb(params,loglikec,lower=rep(1e-10,M+K+L),upper=c(rep(1e10,M+K),rep(1e10,L)),N=N,n=n,M=M,K=K,L=L,mu=mu,MEAN=MEAN,SSE=SSE,pi=pi,classprob=classprob,groupind=groupind)
        a.vec=rep(est$par[1:M],L*K)
        b.vec=rep(rep(est$par[(M+1):(M+K)],each=M),L)
        d.vec=rep(c(est$par[(M+K+1):(M+K+L-1)],1),each=M*K)
        #d.vec=rep(est$par[(M+K+1):(M+K+L)],each=M*K)
        params=est$par
      }
      pi = colSums(classprob)+prior-5/(M*K)
      #pi=colSums(classprob)
      pi = ifelse(pi<0,0,pi); pi=pi/sum(pi);
      mm = post_pi(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE,pi)
      m = mm$pimat
      m.rowsum = rowSums(m)
      loglik[i] = sum(log(m.rowsum))
      classprob = m/m.rowsum
      if(abs(loglik[i]-loglik[i-1])<ltol || loglik[i]<loglik[i-1]) break;      
    }          
  }else{
    if(is.null(pi)){
      pi.a = rep(1,M)/M
      pi.lambda = rep(1,K)/K
      pi.c = c(L,rep(1,L-1))/sum(c(L,rep(1,L-1)))
    }
    loglik = rep(NA,maxiter)
    pi.c = indep_post_pi(N,n,M,K,L,a.vec,b.vec,c.vec,mu,MEAN,SSE, pi.a, pi.lambda, pi.c, prior,group.c)$indep.pi 
    pi.lambda = indep_post_pi(N,n,M,K,L,a.vec,b.vec,c.vec,mu,MEAN,SSE, pi.a, pi.lambda, pi.c, prior,group.lambda)$indep.pi
    m = indep_post_pi(N,n,M,K,L,a.vec,b.vec,c.vec,mu,MEAN,SSE, pi.a, pi.lambda, pi.c,prior, group.a)
    classprob = m$normpi.mat
    pi.a = m$indep.pi
    loglik[1] = sum(log(rowSums(m$pi.mat)))
    pi = rep(pi.a,L*K)*rep(rep(pi.lambda,each=M),L)*rep(pi.c,each=M*K)
    
    for(i in 2:maxiter){
      if(a.lambda.c.est==TRUE){
        #est=nlminb(params,loglike,lower=rep(1e-10,M+K+L),upper=c(rep(100000,M+K),rep(1,L)),N=N,n=n,M=M,K=K,L=L,mu=mu,MEAN=MEAN,SSE=SSE,pi=pi,classprob=classprob,groupind=groupind)
        est=nlminb(params,loglike,gradloglike,lower=rep(1e-10,M+K+L),upper=c(rep(100000,M+K),rep(1,L)),N=N,n=n,M=M,K=K,L=L,mu=mu,MEAN=MEAN,SSE=SSE,pi=pi,classprob=classprob,groupind=groupind)
        a.vec=rep(est$par[1:M],L*K)
        b.vec=rep(rep(est$par[(M+1):(M+K)],each=M),L)
        d.vec=rep(c(est$par[(M+K+1):(M+K+L-1)],1),each=M*K)
        #d.vec=rep(est$par[(M+K+1):(M+K+L)],each=M*K)
        params=est$par
      }
      pi.c = indep_post_pi(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE, pi.a, pi.lambda, pi.c, prior,group.c)$indep.pi
      pi.lambda = indep_post_pi(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE, pi.a, pi.lambda, pi.c,prior, group.lambda)$indep.pi
      m = indep_post_pi(N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE, pi.a, pi.lambda, pi.c,prior, group.a)
      classprob = m$normpi.mat
      pi.a = m$indep.pi
      loglik[i] = sum(log(rowSums(m$pi.mat)))
      pi = rep(pi.a,L*K)*rep(rep(pi.lambda,each=M),L)*rep(pi.c,each=M*K)
      if(abs(loglik[i]-loglik[i-1])<ltol || loglik[i]<loglik[i-1]) break;     
    }   
  }
  lambda.vec=a.vec/b.vec
  if (usePointMass==TRUE){
    c.vec=c((n*d.vec/(1-d.vec))[1:(L-1)],Inf)
  }else{
    c.vec=c((n*d.vec/(1-d.vec))[1:(L-1)],1e10)
  }
  #c.vec=n*d.vec/(1-d.vec)
  converged = (i< maxiter)
  niter= min(c(i,maxiter))
  return(list(pi=pi,pi.a=pi.a,pi.lambda=pi.lambda,pi.c=pi.c,classprob=classprob,loglik.final=loglik,converged=converged,niter=niter,a.vec=a.vec,lambda.vec=lambda.vec,b.vec=b.vec,d.vec=d.vec,c.vec=c.vec))
}

# Compute P(beta|Y,hat(tau)), where hat(tau) is the posterior mean estimate for tau
# ZeroProb=P(beta=0|Y,hat(tau)), PositiveProb=P(beta>0|Y,hat(tau)), NegativeProb=P(beta<0|Y,hat(tau))
CondPostprob = function(pi,tau,gammaa,gammab,gammadens,normmean,normprec,c.vec){
  ZeroProb = rowSums(pi[,c.vec==Inf,drop=FALSE])
  Cond.pi = pi*gammadens/rowSums(pi*gammadens)
  normsd=1/sqrt(outer(tau,normprec))
  PositiveProb = rowSums(Cond.pi[,c.vec==Inf,drop=FALSE]*pnorm(0,mean=normmean[,c.vec==Inf,drop=FALSE],sd=normsd[,c.vec==Inf,drop=FALSE],lower.tail=TRUE))
  NegativeProb =  1-PositiveProb-ZeroProb
  return(list(ZeroProb=ZeroProb,PositiveProb=PositiveProb,NegativeProb=NegativeProb,condpi=Cond.pi,normsd=normsd,normmean=normmean))
}

# Compute P(beta|Y)
# ZeroProb=P(beta=0|Y), PositiveProb=P(beta>0|Y), NegativeProb=P(beta<0)
Postprob = function(mu,pi,gammaa,gammab,normmean,normc,c.vec){
  ZeroProb = rowSums(pi[,c.vec==Inf,drop=FALSE])
  mumat=outer(mu,rep(1,length(c.vec)))
  T.std = ((mumat-normmean)/sqrt(gammab/(gammaa*normc)))[,c.vec!=Inf,drop=FALSE]
  pval.mat=2*pt(abs(T.std),df=gammaa[,c.vec!=Inf]*2,lower.tail=FALSE)
  #pval.mat=rowSums(pi*pval.mat)
  PositiveProb = rowSums(pi[,c.vec!=Inf,drop=FALSE]*pt(T.std,df=gammaa[,c.vec!=Inf,drop=FALSE]*2,lower.tail=TRUE))
  NegativeProb =  1-PositiveProb-ZeroProb
  return(list(ZeroProb=ZeroProb,PositiveProb=PositiveProb,NegativeProb=NegativeProb,pval.mat=pval.mat))
}

# Compute local FALSE sign rate & FALSE discovery rate
computefdr = function(ZeroProb,PositiveProb,NegativeProb){
  localfsr = ifelse(PositiveProb<NegativeProb,PositiveProb+ZeroProb,NegativeProb+ZeroProb)
  if(sum(ZeroProb)>0){
    localfdr = ZeroProb
  }else{
    localfdr = ifelse(PositiveProb < NegativeProb, 2 * PositiveProb + ZeroProb, 
                      2 * NegativeProb + ZeroProb)
  }  
  return(list(localfsr=localfsr,localfdr=localfdr))
}

lmreg=function(Y,fac,nA,nB){
  nfac=length(unique(fac))
  design=model.matrix(~fac)
  N=dim(Y)[1]
  betahat=rep(NA,N)
  sebetahat=rep(NA,N)
  for (i in 1:N){
    resp = Y[i,]
    fit = lm(resp~0+design)
    betahat[i]=fit$coef[2]
    sebetahat[i]=coef(summary(fit))[, "Std. Error"][2]
  }
  return(list(betahat=betahat,sebetahat=sebetahat))
}

lmreg1=function(Y,fac,nA,nB){
  betahat=apply(Y[,fac==2],1,mean)-apply(Y[,fac==1],1,mean)
  Yfit=cbind(outer(apply(Y[,fac==1],1,mean),rep(1,nA)),outer(apply(Y[,fac==2],1,mean),rep(1,nA)))
  design=model.matrix(~fac)
  XXinv=solve(t(design)%*%design)[2,2]
  #SSE=sigmahat2*(nA+nB-2)/XXinv
  SSE=apply((Y-Yfit)^2,1,sum)
  MEAN=betahat/sqrt(XXinv*(nA+nB-1))
  return(list(MEAN=MEAN,SSE=SSE,scale=sqrt(XXinv*(nA+nB-1))))
}
###########NEED TESTING!#########
autoselect.precMulti = function(Y){
  n = dim(Y)[2]
  Y.var = apply(Y,1,var)*(n-1)/n
  lmin = 1/median(Y.var)/4
  lmax = lmin*64
  npoint = ceiling(log2(lmax/lmin))
  precMulti = 2^((-npoint):0) * lmax
  return(precMulti) 
}

# Sample posterior beta & tau from P(beta|Y), P(tau|Y)
# obs: obs need to be sampled
# nsamp: num of sampled beta & tau for each obs
posterior_sample_jash = function(post,nsamp,obs){
  component.tau = as.vector(apply(post$pi[obs,],1,sample_component,nsamp=nsamp))
  ML = ncol(post$pi)
  obs.vec = rep(obs,each=nsamp)
  index.tau = cbind(obs.vec,component.tau)
  tau = rgamma(length(component.tau),shape=post$gammaa[index.tau],rate=post$gammab[index.tau])
  normsd=1/sqrt(outer(tau,post$normprec))
  beta = rnorm(length(tau),mean=post$normmean[index.tau],sd=normsd[cbind(1:length(tau),component.tau)])
  res = list(beta=matrix(beta,ncol=nsamp,byrow=TRUE),tau=matrix(tau,ncol=nsamp,byrow=TRUE),obs.vec=obs.vec)
  return(res)
}

# Y: data matrix, N by n
# precShape: vector of a_m
# precMulti: vector of lambda_k
# compprecPrior: vector of c_l
# fac: group factor vec
jash = function(Y, fac, auto=FALSE, precShape=NULL, precMulti=NULL, compprecPrior=NULL, mu=NULL, pi=NULL,prior=NULL, usePointMass=FALSE, localfdr=TRUE, a.lambda.c.est=TRUE, SGD=FALSE){
  N = dim(Y)[1]
  
  if(is.null(mu)){
    mu = rep(0,N)
  }
  if(auto==TRUE){
    precMulti = autoselect.precPrior(Y)
  }
  if(is.null(precShape)){
    precShape = c(0.01,0.1,1,10,100)
  }
  if(is.null(precMulti)){
    precMulti = c(0.25,0.5,1,2,4)
  }
  if(is.null(compprecPrior) && usePointMass==FALSE){
    compprecPrior = c(1,1e10)
  }else if(is.null(compprecPrior) && usePointMass==TRUE){
    compprecPrior=c(1,Inf)
  }else if(usePointMass==TRUE){
    compprecPrior=c(compprecPrior,Inf)
  }
  if(a.lambda.c.est==TRUE){
    precMulti=1
    precShape=1
  }
  
  if(length(unique(fac))==1){
    n = dim(Y)[2]
    MEAN = apply(Y,1,mean)
    SSE = apply(Y,1,var)*(n-1) 
    reg=list(scale=1)
  }else{
    nA=sum(fac==1)
    nB=sum(fac==2)
    reg=lmreg1(Y,fac,nA,nB)
    MEAN=reg$MEAN
    SSE=reg$SSE
    n=nA+nB-1
  }
  
  M = length(precShape)
  K = length(precMulti)
  L = length(compprecPrior) 
  params=c(precShape,precShape/precMulti,(compprecPrior/(n+compprecPrior))[1:(L-1)])
  a.vec = rep(precShape,L*K)
  lambda.vec = rep(rep(precMulti,each=M),L)
  b.vec = a.vec/lambda.vec
  c.vec = rep(compprecPrior,each=M*K)
  d.vec=1/(n/c.vec+1)
  group.a = rep(1:M,L*K)
  group.lambda = rep(rep(1:K,each=M),L)
  group.c = rep(1:L,each=M*K)
    
  pifit = EMest_pi(params,N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE, pi,prior,a.lambda.c.est,SGD,indepprior=FALSE,ltol=0.0001, maxiter=2000, usePointMass)
  post = post_distn(N,n,M,K,L,pifit$a.vec,pifit$b.vec,pifit$c.vec,mu,MEAN,SSE,pifit$pi)
  #condpost = CondPostprob(post$pi,post$tau,post$gammaa,post$gammab,post$gammadens,post$normmean,post$normprec,post$c.vec)
  postprob = Postprob(mu,post$pi,post$gammaa,post$gammab,post$normmean,post$normc,pifit$c.vec)
  if(localfdr==TRUE){
    localfdr = computefdr(postprob$ZeroProb,postprob$PositiveProb,postprob$NegativeProb)$localfdr
    qvalue = qval.from.localfdr(localfdr)
  }else{
    localfdr = NULL
    qvalue = NULL
  }
  if(sum(post$c.vec==Inf)){
    null.postprob = postprob$ZeroProb
  }else{
    null.postprob = NULL
  }
  #fitted = list(pi=pifit$pi,a.vec=post$a.vec,lambda.vec=post$lambda.vec,c.vec=post$c.vec,pi.a=pifit$pi.a,pi.lambda=pifit$pi.lambda,pi.c=pifit$pi.c)
  PosteriorMean=post$beta*reg$scale
  return(list(PosteriorMean=PosteriorMean,PosteriorPrec=post$tau,pifit=pifit,post=post,postprob=postprob,localfdr=localfdr,qvalue=qvalue,null.postprob=null.postprob,
              a.vec=pifit$a.vec,lambda.vec=pifit$lambda.vec,c.vec=pifit$c.vec,mu=mu))
}

jasha=function(betahat,betahatsd,df,auto=FALSE, precShape=NULL, precMulti=NULL, compprecPrior=NULL, mu=NULL, pi=NULL,prior=NULL, usePointMass=FALSE, localfdr=TRUE, a.lambda.c.est=TRUE, SGD=FALSE){
  N = length(betahat)
  n=df+1
  
  if(is.null(mu)){
    mu = rep(0,N)
  }
  if(auto==TRUE){
    precMulti = autoselect.precPrior(Y)
  }
  if(is.null(precShape)){
    precShape = c(0.01,0.1,1,10,100)
  }
  if(is.null(precMulti)){
    precMulti = c(0.25,0.5,1,2,4)
  }
  if(is.null(compprecPrior) && usePointMass==FALSE){
    compprecPrior = c(1,1e10)
  }else if(is.null(compprecPrior) && usePointMass==TRUE){
    compprecPrior=c(1,Inf)
  }else if(usePointMass==TRUE){
    compprecPrior=c(compprecPrior,Inf)
  }
  if(a.lambda.c.est==TRUE){
    precMulti=1
    precShape=1
  }
  
  MEAN=betahat/sqrt(n)
  SSE=betahatsd^2*n
  
  M = length(precShape)
  K = length(precMulti)
  L = length(compprecPrior) 
  params=c(precShape,precShape/precMulti,(compprecPrior/(n+compprecPrior))[1:(L-1)])
  a.vec = rep(precShape,L*K)
  lambda.vec = rep(rep(precMulti,each=M),L)
  b.vec = a.vec/lambda.vec
  c.vec = rep(compprecPrior,each=M*K)
  d.vec=1/(n/c.vec+1)
  group.a = rep(1:M,L*K)
  group.lambda = rep(rep(1:K,each=M),L)
  group.c = rep(1:L,each=M*K)
  
  pifit = EMest_pi(params,N,n,M,K,L,a.vec,b.vec,d.vec,mu,MEAN,SSE, pi,prior,a.lambda.c.est,SGD,indepprior=FALSE,ltol=0.0001, maxiter=2000, usePointMass)
  post = post_distn(N,n,M,K,L,pifit$a.vec,pifit$b.vec,pifit$c.vec,mu,MEAN,SSE,pifit$pi)
  #condpost = CondPostprob(post$pi,post$tau,post$gammaa,post$gammab,post$gammadens,post$normmean,post$normprec,post$c.vec)
  postprob = Postprob(mu,post$pi,post$gammaa,post$gammab,post$normmean,post$normc,pifit$c.vec)
  if(localfdr==TRUE){
    localfdr = computefdr(postprob$ZeroProb,postprob$PositiveProb,postprob$NegativeProb)$localfdr
    qvalue = qval.from.localfdr(localfdr)
  }else{
    localfdr = NULL
    qvalue = NULL
  }
  if(sum(post$c.vec==Inf)){
    null.postprob = postprob$ZeroProb
  }else{
    null.postprob = NULL
  }
  
  #fitted = list(pi=pifit$pi,a.vec=post$a.vec,lambda.vec=post$lambda.vec,c.vec=post$c.vec,pi.a=pifit$pi.a,pi.lambda=pifit$pi.lambda,pi.c=pifit$pi.c)
  PosteriorMean=post$beta*sqrt(n)
  params=c(pifit$a.vec[1],pifit$b.vec[1]/n,pifit$c.vec[1])
  loglik=loglike(params,N,n,M,K,L,mu,MEAN,rep(0,N),pifit$pi,0)
  
  return(list(PosteriorMean=PosteriorMean,PosteriorPrec=post$tau,pifit=pifit,post=post,postprob=postprob,localfdr=localfdr,qvalue=qvalue,null.postprob=null.postprob,
              a.vec=pifit$a.vec,lambda.vec=pifit$lambda.vec,c.vec=pifit$c.vec,mu=mu,loglik=loglik))
}