# Y: N by n observation matrix (N tests, n replicate obs for each test)
# Model: Y_{jk}|beta_j,tau_j ~ N(beta_j,1/tau_j)
# Prior: beta_j|tau_j~N(mu_j,1/(c_l*tau_j)), tau_j~lambda_k*Gamma(a_m,a_m), w.p. pi_{klm}
# Mixture constraint: sum(pi_{klm})=1
# Or: sum(pi_k)=1, sum(pi_l)=1, sum(pi_m)=1 for each param
# a_m, lambda_k, c_l, mu_j are known
######################################################

# post_pi computes the posterior probability 
# P(beta_j,tau_j belongs to the i'th component | Y)
post_pi = function(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst,pi){
  post.gammab=0.5*(outer(Y.sst,rep(1,M*L*K))+outer((Y.mean-mu)^2,n/(n/c.vec+1))+outer(rep(1,N),2*a.vec/lambda.vec))
  post.pi.mat = exp(outer(rep(1,N),log(pi)+a.vec*log(a.vec/lambda.vec)-lgamma(a.vec)-0.5*log(n/c.vec+1)+lgamma(a.vec+n/2))-
                      log(post.gammab)*
                      outer(rep(1,N),a.vec+n/2))
  return(list(gammab=post.gammab,pimat=post.pi.mat))
  
#   post.pi.mat = exp(outer(rep(1,N),log(pi)+a.vec*log(a.vec/lambda.vec)-lgamma(a.vec)+log(sqrt(1/(n/c.vec+1)))+lgamma(a.vec+n/2))-
#     log(0.5*(outer(Y.sst,rep(1,M*L*K))+outer((Y.mean-mu)^2,n/(n/c.vec+1))+outer(rep(1,N),2*a.vec/lambda.vec)))*
#     outer(rep(1,N),a.vec+n/2))
 
#   post.pi.mat = outer(rep(1,N),pi*(a.vec/lambda.vec)^a.vec/gamma(a.vec)*sqrt(1/(n/c.vec+1))*(gamma(a.vec+n/2)))/
#     (0.5*(outer(Y.sst,rep(1,M*L*K))+outer((Y.mean-mu)^2,n/(n/c.vec+1))+outer(rep(1,N),2*a.vec/lambda.vec)))^
#     outer(rep(1,N),a.vec+n/2)
  return(post.pi.mat)
}

# indep_post_pi computes the posterior probability 
# P(beta_j,tau_j belongs to the i'th component of precShape/precMulti.compprecPrior| Y)
# groupind: i indicates that the current component correspond to the i'th comp of precShape/precMulti/compprecPrior
indep_post_pi = function(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi.a, pi.lambda, pi.c, prior, groupind){
  pi = rep(pi.a,L*K)*rep(rep(pi.lambda,each=M),L)*rep(pi.c,each=M*K)
  mm=post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst,pi)
  pi.mat = mm$pimat
  postb=mm$gammab
  norm.pi.mat = normalized_pi(pi.mat)
  indep.post.pi = rep(0,max(groupind))
  temp = colSums(norm.pi.mat)+prior-5/(M*K)
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
post_distn = function(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst,pi){ 
  post.pi = normalized_pi(post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst,pi)$pimat)

  post.gamma.parama = outer(rep(1,N),a.vec+n/2)
  post.gamma.paramb = 0.5*(outer(Y.sst,rep(1,M*L*K))+outer((Y.mean-mu)^2,n/(n/c.vec+1))+outer(rep(1,N),2*a.vec/lambda.vec))
  
  # Posterior mean: E(tau|Y)
  post.tau = apply(post.gamma.parama/post.gamma.paramb*post.pi,1,sum)/apply(post.pi,1,sum)
  
  # Posterior mean: E(beta|Y)
  post.norm.mean = outer(Y.mean,n/(n+c.vec))+outer(mu,1/(n/c.vec+1))
  post.norm.prec = n+c.vec # if tau given, outer(tau,n+c.vec)
  post.norm.c = outer(rep(1,N),n+c.vec)
  post.gamma.dens = dgamma(outer(post.tau,rep(1,M*L*K)),shape=post.gamma.parama,rate=post.gamma.paramb)
  post.beta = apply(post.pi*post.gamma.dens*post.norm.mean,1,sum)/
    apply(post.pi*post.gamma.dens,1,sum)
  
  return(list(pi=post.pi,tau=post.tau,beta=post.beta,gammaa=post.gamma.parama,gammab=post.gamma.paramb,
              gammadens=post.gamma.dens,normmean=post.norm.mean,normc=post.norm.c,normprec=post.norm.prec,c.vec=c.vec,a.vec=a.vec,lambda.vec=lambda.vec))
}

# Use EM algorithm to estimate a, lambda
# lik: likelihood matrix
# SGD: use stochastic gradient descent method
a_lambda_est=function(classprob,postb,N,n,M,K,L,a.vec,lambda.vec,c.vec,group.a,group.lambda,SGD,learnrate){
  b.vec=a.vec/lambda.vec

  if(SGD==TRUE){    
    dl_da=sum(classprob*(outer(rep(1,N),log(b.vec))-log(postb)+outer(rep(1,N),digamma(a.vec+n/2)-digamma(a.vec))))
    dl_db=sum(classprob*(outer(rep(1,N),a.vec/b.vec)-outer(rep(1,N),a.vec+n/2)/postb))
    a.vec=pmin(a.vec-learnrate*dl_da,1e-10)
    b.vec=pmin(b.vec-learnrate*dl_db,1e-10) 
    print(dl_da)
    print(dl_db)
  }else{
    for(i in 1:max(group.a)){
      idx=(group.a==i)
      fa=function(a) {        
        res=sum(classprob[idx,]*(outer(rep(1,N),log(b.vec[idx]))-log(postb[idx])+outer(rep(1,N),rep(digamma(a+n/2)-digamma(a),length(idx)))))
        return(res)
      }
      a.vec[idx]=rep(uniroot(fa,c(1e-10,100000))$root,length(idx))    
    }
    for(i in 1:max(group.lambda)){
      idx=(group.lambda==i)
      fb=function(b) {
        res=sum(classprob[idx]*(outer(rep(1,N),a.vec[idx]/b)-outer(rep(1,N),a.vec[idx]+n/2)/postb[idx]))
        return(res)
      }
      b.vec[idx]=rep(uniroot(fb,c(1e-10,100000))$root,length(idx))    
    }
  }
  lambda.vec=a.vec/b.vec
  return(list(a.vec=a.vec,lambda.vec=lambda.vec))
}



# Use EM algorithm to estimate pi from data
# ltol: likelihood convergence tolerance
# maxiter: max number of iterations
# indepprior: whether assume that there are independent priors for a_m, lambda_k and c_l
# i.e. pi_{klm}=pi.lambda_k*pi.c_l*pi.a_m, where sum(pi.a_m)=sum(pi.lambda_k)=sum(pi.c_l)
EMest_pi = function(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi, prior, a.lambda.est, SGD,indepprior=TRUE,ltol=0.0001, maxiter=2000){
  group.a = rep(1:M,L*K)
  group.lambda = rep(rep(1:K,each=M),L)
  group.c = rep(1:L,each=M*K)
  
  if(is.null(prior)){ # NEED TESTING
    prior=rep(5/(L*M*K),L*M*K)
    prior[c.vec==max(c.vec)] = 5/(M*K)     
  }else if(prior=="uniform"){
    prior=rep(5/(M*K),M*L*K)
  }
  if(indepprior==FALSE){
    pi.a=NULL; pi.lambda=NULL; pi.c=NULL
    if(is.null(pi)){
      pi=rep(1,K*L*M)/(K*L*M)
      pi[c.vec==max(c.vec)]=L
      pi=pi/sum(pi)
    }
    loglik = rep(NA,maxiter)
    
    mm = post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst,pi)
    m = mm$pimat
    m.rowsum = rowSums(m)
    loglik[1] = sum(log(m.rowsum))
    classprob = m/m.rowsum
    
    for(i in 2:maxiter){
      if(a.lambda.est==TRUE){
        est=a_lambda_est(classprob,mm$gammab,N,n,M,K,L,a.vec,lambda.vec,c.vec,group.a,group.lambda,SGD,learnrate=1/(100*i))
        a.vec=est$a.vec
        lambda.vec=est$lambda.vec
      }
      pi = colSums(classprob)+prior-5/(M*K)
      pi = ifelse(pi<0,0,pi); pi=pi/sum(pi);
      mm = post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst,pi)
      m = mm$pimat
      m.rowsum = rowSums(m)
      loglik[i] = sum(log(m.rowsum))
      classprob = m/m.rowsum
      if(abs(loglik[i]-loglik[i-1])<ltol) break;      
    }          
  }else{
    if(is.null(pi)){
      pi.a = rep(1,M)/M
      pi.lambda = rep(1,K)/K
      pi.c = c(L,rep(1,L-1))/sum(c(L,rep(1,L-1)))
    }
    loglik = rep(NA,maxiter)
    pi.c = indep_post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi.a, pi.lambda, pi.c, prior,group.c)$indep.pi 
    pi.lambda = indep_post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi.a, pi.lambda, pi.c, prior,group.lambda)$indep.pi
    m = indep_post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi.a, pi.lambda, pi.c,prior, group.a)
    classprob = m$normpi.mat
    pi.a = m$indep.pi
    loglik[1] = sum(log(rowSums(m$pi.mat)))
    
    for(i in 2:maxiter){
      if(a.lambda.est==TRUE){
        est=a_lambda_est(classprob,m$postb,N,n,M,K,L,a.vec,lambda.vec,c.vec,group.a,group.lambda,SGD,learnrate=1/(i^2))
        a.vec=est$a.vec
        lambda.vec=est$lambda.vec
      }
      pi.c = indep_post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi.a, pi.lambda, pi.c, prior,group.c)$indep.pi
      pi.lambda = indep_post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi.a, pi.lambda, pi.c,prior, group.lambda)$indep.pi
      m = indep_post_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi.a, pi.lambda, pi.c,prior, group.a)
      classprob = m$normpi.mat
      pi.a = m$indep.pi
      loglik[i] = sum(log(rowSums(m$pi.mat)))
      if(abs(loglik[i]-loglik[i-1])<ltol) break;     
    }   
    pi = rep(pi.a,L*K)*rep(rep(pi.lambda,each=M),L)*rep(pi.c,each=M*K)
  }

  converged = (i< maxiter)
  niter= min(c(i,maxiter))
  return(list(pi=pi,pi.a=pi.a,pi.lambda=pi.lambda,pi.c=pi.c,classprob=classprob,loglik.final=loglik,converged=converged,niter=niter,a.vec=a.vec,lambda.vec=lambda.vec))
}

# Compute P(beta|Y,hat(tau)), where hat(tau) is the posterior mean estimate for tau
# ZeroProb=P(beta=0|Y,hat(tau)), PositiveProb=P(beta>0|Y,hat(tau)), NegativeProb=P(beta<0|Y,hat(tau))
CondPostprob = function(pi,tau,gammaa,gammab,gammadens,normmean,normprec,c.vec){
  ZeroProb = rowSums(pi[,c.vec==Inf,drop=FALSE])
  Cond.pi = pi*gammadens/rowSums(pi*gammadens)
  normsd=1/sqrt(outer(tau,normprec))
  PositiveProb = rowSums(Cond.pi*pnorm(0,mean=normmean,sd=normsd,lower.tail=TRUE))
  NegativeProb =  1-PositiveProb-ZeroProb
  return(list(ZeroProb=ZeroProb,PositiveProb=PositiveProb,NegativeProb=NegativeProb,condpi=Cond.pi,normsd=normsd,normmean=normmean))
}

# Compute P(beta|Y)
# ZeroProb=P(beta=0|Y), PositiveProb=P(beta>0|Y), NegativeProb=P(beta<0)
Postprob = function(pi,gammaa,gammab,normmean,normc,c.vec){
  ZeroProb = rowSums(pi[,c.vec==Inf,drop=FALSE])
  T.std = (-normmean/sqrt(gammab/(gammaa*normc)))[,c.vec!=Inf]
  PositiveProb = rowSums(pi[,c.vec!=Inf]*pt(T.std,df=gammaa[,c.vec!=Inf]*2,lower.tail=TRUE))
  NegativeProb =  1-PositiveProb-ZeroProb
  return(list(ZeroProb=ZeroProb,PositiveProb=PositiveProb,NegativeProb=NegativeProb))
}

# Compute local false sign rate & false discovery rate
computefdr = function(ZeroProb,PositiveProb,NegativeProb){
  localfsr = ifelse(PositiveProb<NegativeProb,PositiveProb+ZeroProb,NegativeProb+ZeroProb)
  if(sum(ZeroProb)>0){
    localfdr = localfsr
  }else{
    localfdr = 2*localfsr
  }  
  return(list(localfsr=localfsr,localfdr=localfdr))
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
jash = function(Y, auto=FALSE, precShape=NULL, precMulti=NULL, compprecPrior=NULL, mu=NULL, pi=NULL,prior=NULL, usePointMass=FALSE, localfdr=TRUE, a.lambda.est=TRUE, SGD=FALSE){
  N = dim(Y)[1]
  n = dim(Y)[2]
  
  if(is.null(mu)){
    mu = rep(0,N)
  }
  if(auto==TRUE){
    precMulti = autoselect.precMulti(Y)
  }
  if(is.null(precShape)){
    precShape = c(0.01,0.1,1,10,100)
  }
  if(is.null(precMulti)){
    precMulti = c(0.25,0.5,1,2,4)
  }
  if(is.null(compprecPrior) && usePointMass==FALSE){
    compprecPrior = c(1,1000000)
  }else if(is.null(compprecPrior) && usePointMass==TRUE){
    compprecPrior=c(1,Inf)
  }else if(usePointMass==TRUE){
    compprecPrior=c(compprecPrior,Inf)
  }
  if(a.lambda.est==TRUE){
    precMulti=1
    precShape=1
  }
  
  M = length(precShape)
  K = length(precMulti)
  L = length(compprecPrior) 
  a.vec = rep(precShape,L*K)
  lambda.vec = rep(rep(precMulti,each=M),L)
  c.vec = rep(compprecPrior,each=M*K)
  group.a = rep(1:M,L*K)
  group.lambda = rep(rep(1:K,each=M),L)
  group.c = rep(1:L,each=M*K)
  Y.mean = apply(Y,1,mean)
  Y.sst = apply(Y,1,var)*(n-1)
  
  pifit = EMest_pi(N,n,M,K,L,a.vec,lambda.vec,c.vec,mu,Y.mean,Y.sst, pi,prior,a.lambda.est,SGD,indepprior=TRUE,ltol=0.0001, maxiter=2000)
  post = post_distn(N,n,M,K,L,pifit$a.vec,pifit$lambda.vec,c.vec,mu,Y.mean,Y.sst,pifit$pi)
  #condpost = CondPostprob(post$pi,post$tau,post$gammaa,post$gammab,post$gammadens,post$normmean,post$normprec,post$c.vec)
  postprob = Postprob(post$pi,post$gammaa,post$gammab,post$normmean,post$normc,post$c.vec)
  if(localfdr==TRUE){
    localfdr = computefdr(postprob$ZeroProb,postprob$PositiveProb,postprob$NegativeProb)$localfdr
    qvalue = qval.from.lfdr(localfdr)
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
  
  return(list(PosteriorMean=post$beta,PosteriorPrec=post$tau,pifit=pifit,post=post,postprob=postprob,localfdr=localfdr,qvalue=qvalue,null.postprob=null.postprob,
              a.vec=post$a.vec,lambda.vec=post$lambda.vec,c.vec=post$c.vec,mu=mu))
}

