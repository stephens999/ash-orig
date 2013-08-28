source('ash.R')

# Y: N by n observation matrix (N tests, n replicate obs for each test)
# Model: Y_{jk}|beta_j,tau_j ~ N(beta_j,1/tau_j)
# Prior: beta_j|tau_j~N(mu_j,1/(c_l*tau_j)), tau_j~Gamma(a_m,b_m), w.p. pi_{ml}
# Mixture constraint: sum(pi_{ml})=1
# a_m, b_m, c_l, mu_j are known
######################################################

# post_pi computes the posterior probability 
# P(beta_j,tau_j belongs to the i'th component | Y)
post_pi = function(Y, precPrior, compprecPrior, mu, pi){
  N = dim(Y)[1]
  n = dim(Y)[2]
  M = dim(precPrior)[1]
  L = length(compprecPrior)
  
  a.vec = rep(precPrior[,1],L)
  b.vec = rep(precPrior[,2],L)
  c.vec = rep(compprecPrior,each=M)
  Y.mean = apply(Y,1,mean)
  Y.sst = apply(Y,1,var)*(n-1)
  
  post.pi.mat = outer(rep(1,N),pi*b.vec^a.vec/gamma(a.vec)*sqrt(1/(n/c.vec+1))*gamma(a.vec+n/2))/
    (0.5*(outer(Y.sst,rep(1,M*L))+outer((Y.mean-mu)^2,n/(n/c.vec+1))+outer(rep(1,N),2*b.vec)))^
    outer(rep(1,N),a.vec+n/2)
  return(post.pi.mat)
}

# Normalize pi to make sum(pi)=1
normalized_pi = function(pi.mat){
  n = dim(pi.mat)[2]
  pi.normalized = pi.mat/outer(rowSums(pi.mat),rep(1,n))  
  return(pi.normalized)
}

# Compute posterior distribution P(tau|Y), P(beta|Y,tau)
# precPrior: M by 2 matrix, all pairs of (a_m, b_m) (a_m>0, b_m>0)
# compprecPrior: L vector, all c_l (c_l>=0)
# pi: (M*L) vector
# mu: N vector
post_distn = function(Y, precPrior, compprecPrior, mu, pi){
  N = dim(Y)[1]
  n = dim(Y)[2]
  M = dim(precPrior)[1]
  L = length(compprecPrior)
  
  a.vec = rep(precPrior[,1],L)
  b.vec = rep(precPrior[,2],L)
  c.vec = rep(compprecPrior,each=M)
  Y.mean = apply(Y,1,mean)
  Y.sst = apply(Y,1,var)*(n-1)
  
  post.pi = post_pi(Y, precPrior, compprecPrior, mu, pi)

  post.gamma.parama = outer(rep(1,N),a.vec+n/2)
  post.gamma.paramb = 0.5*(outer(Y.sst,rep(1,M*L))+outer((Y.mean-mu)^2,n/(n/c.vec+1))+outer(rep(1,N),2*b.vec))
  
  # Posterior mean: E(tau|Y)
  post.tau = apply(post.gamma.parama/post.gamma.paramb*post.pi,1,sum)/apply(post.pi,1,sum)
  
  # Posterior mean: E(beta|Y)
  post.norm.mean = outer(Y.mean,n/(n+c.vec))+outer(mu,1/(n/c.vec+1))
  post.norm.prec = n+c.vec # if tau given, outer(tau,n+c.vec)     
  post.gamma.dens = dgamma(outer(post.tau,rep(1,M*L)),shape=post.gamma.parama,rate=post.gamma.paramb)
  post.beta = apply(post.pi*post.gamma.dens*post.norm.mean,1,sum)/
    apply(post.pi*post.gamma.dens,1,sum)
  
  return(list(pi=post.pi,tau=post.tau,beta=post.beta,gammaa=post.gamma.parama,gammab=post.gamma.paramb,
              gammadens=post.gamma.dens,normmean=post.norm.mean,normprec=post.norm.prec,c.vec=c.vec,a.vec=a.vec,b.vec=b.vec))
}

# Use EM algorithm to estimate pi from data
# ltol: likelihood convergence tolerance
# maxiter: max number of iterations
EMest_pi = function(Y, precPrior, compprecPrior, mu, prior, ltol=0.0001, maxiter=2000){
  varvec = c(outer(1+1/compprecPrior,1/(precPrior[,1]*precPrior[,2]))) # Var(y_{jk}) for each comp (when a_m>1)
  null.comp = which.min(varvec)
  ML = length(varvec)
  if(is.null(prior)){
    prior=rep(1/(ML-1),ML)
    prior[null.comp] = 1
  }else if(prior=="uniform"){
    prior = rep(1,ML)
  }
  prior=prior/sum(prior)
  
  loglik = rep(NA,maxiter)
  loglik[1] = -Inf
  pi = apply(normalized_pi(post_pi(Y, precPrior, compprecPrior, mu, prior)),2,mean)
  for(i in 2:maxiter){
    m = post_pi(Y, precPrior, compprecPrior, mu, pi)
    m.rowsum = rowSums(m)
    loglik[i] = sum(log(m.rowsum))
    classprob = m/m.rowsum
    pi = apply(classprob,2,mean)
    if(abs(loglik[i]-loglik[i-1])<ltol) break;
    
  }
  converged = (i< maxiter)
  niter= min(c(i,maxiter))
  return(list(pi=pi,classprob=classprob,loglik.final=loglik,converged=converged,niter=niter,varvec=varvec))
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

# Compute local false sign rate & false discovery rate
computefdr = function(ZeroProb,PositiveProb,NegativeProb){
  localfsr = ifelse(PositiveProb<NegativeProb,PositiveProb,NegativeProb)+ZeroProb
  if(sum(ZeroProb)>0){
    localfdr = localfsr
  }else{
    localfdr = 2*localfsr
  }  
  return(list(localfsr=localfsr,localfdr=localfdr))
}

###########NEED TESTING!#########
autoselect.precPrior = function(Y){
  Y.mean = apply(Y,1,mean)
  Y.sd = apply(Y,1,sd)
  bmin = 1/max(Y.sd)
  bmax = bmin*10
  npoint = ceiling(log2(bmax/bmin))
  b.vec = 2^((-npoint):0) * bmax
  #bmin = 1/max(Y.sd)/10
  #bmax = bmin*100
  #npoint = ceiling(log(bmax/bmin)/log(3))
  #b.vec = 3^((-npoint):0) * bmax
  a.vec = rep(2,length(b.vec))
  precPrior = matrix(c(a.vec,b.vec),ncol=2)
  return(precPrior) 
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
#   gammadens = dgamma(outer(tau,rep(1,ML)),shape=post$gammaa[obs.vec,],rate=post$gammab[obs.vec,])
#   Cond.pi = post$pi[obs.vec,]*gammadens/rowSums(post$pi[obs.vec,]*gammadens)
#   component.beta = as.vector(apply(Cond.pi,1,sample_component,nsamp=1))
#   index.beta = cbind(obs.vec,component.tau)
#   normsd=1/sqrt(outer(tau,post$normprec))
#   beta = rnorm(length(component.beta),mean=post$normmean[index.beta],sd=normsd[cbind(1:length(component.beta),component.tau)])
  res = list(beta=matrix(beta,ncol=nsamp,byrow=TRUE),tau=matrix(tau,ncol=nsamp,byrow=TRUE),obs.vec=obs.vec)
  return(res)
}

jash = function(Y, auto=FALSE, precPrior=NULL, compprecPrior=NULL, mu=NULL, prior=NULL, usePointMass=FALSE, localfdr=TRUE){
  N = dim(Y)[1]
  n = dim(Y)[2]
  
  if(is.null(mu)){
    mu = rep(0,N)
  }
  if(auto==TRUE){
    precPrior = autoselect.precPrior(Y)
  }
  if(is.null(precPrior)){
    precPrior = matrix(c(seq(5,2),c(50,20,1,0.1)),ncol=2)
    #precPrior = matrix(c(seq(5,2),c(500,20,1,0.1)),ncol=2)
  }
  if(is.null(compprecPrior)){
    compprecPrior = c(10,1,0.1,0.01)
  }    
  if(usePointMass==TRUE){
    compprecPrior=c(Inf,compprecPrior)
  }
  
  pifit = EMest_pi(Y, precPrior, compprecPrior, mu, prior)
  post = post_distn(Y, precPrior, compprecPrior, mu, pifit$pi)
  condpost = CondPostprob(post$pi,post$tau,post$gammaa,post$gammab,post$gammadens,post$normmean,post$normprec,post$c.vec)
  if(localfdr==TRUE){
    localfdr = computefdr(condpost$ZeroProb,condpost$PositiveProb,condpost$NegativeProb)$localfdr
    qvalue = qval.from.localfdr(localfdr)
  }else{
    localfdr = NULL
    qvalue = NULL
  }
  if(sum(post$c.vec==Inf)){
    null.postprob = condpost$ZeroProb
  }else{
    null.postprob = NULL
  }
  fitted = list(pi=pifit$pi,a.vec=post$a.vec,b.vec=post$b.vec,c.vec=post$c.vec,sdvec=sqrt(pifit$varvec))
  
  return(list(PosteriorMean=post$beta,fitted=fitted,pifit=pifit,post=post,condpost=condpost,localfdr=localfdr,qvalue=qvalue,null.postprob=null.postprob,
              a.vec=post$a.vec,b.vec=post$b.vec,c.vec=post$c.vec,mu=mu))
}

