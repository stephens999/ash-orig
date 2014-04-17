qval.from.lfdr = function(lfdr){
  o = order(lfdr)
  qvalue=rep(NA,length(lfdr))
  qvalue[o] = (cumsum(sort(lfdr))/(1:sum(!is.na(lfdr))))
  return(qvalue)
}

normalize = function(x){return(x/sum(x))}

# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
# mult is the multiplier by which the sds differ across the grid
autoselect.mixsd = function(betahat,sebetahat,mult){
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  } else {
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}

# prior.se of se: se|pi,alpha.vec,c ~ pi*IG(alpha_k,c*(alpha_k-1))
# Likelihood: varhat|se ~ sj*Gamma(n/2,n/2)
# pi, alpha.vec, c: known
# Posterior weight of P(se|varhat) (IG mixture distn)
post_pi = function(n,varhat,alpha.vec,modalpha.vec,c,pi){
  N = length(varhat)
  K = length(alpha.vec)
  post.pi.mat = exp(outer(rep(1,N),log(pi))+n/2*log(n/2)-lgamma(n/2)
                  +(n/2-1)*outer(log(varhat),rep(1,K))
                  +outer(rep(1,N),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+n/2))
                  -outer(rep(1,N),alpha.vec+n/2)*log(outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  return(pimat=post.pi.mat)
}

# Normalize pi to make sum(pi)=1
normalized_pi = function(pi.mat){
  n = dim(pi.mat)[2]
  pi.normalized = pi.mat/outer(rowSums(pi.mat),rep(1,n))  
  return(pi.normalized)
}

# Posterior distn P(se|varhat)
# alpha.vec, c, pi: known
# varhat: observed (standard errors)^2
post_distn_se = function(n,varhat,alpha.vec,modalpha.vec,c,pi){ 
  N = length(varhat)
  K = length(alpha.vec)
  post.pi = normalized_pi(post_pi(n,varhat,alpha.vec,modalpha.vec,c,pi))
  
  post.IG.parama = outer(rep(1,N),alpha.vec+n/2)
  post.IG.paramb = outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))
  
  # Posterior mean: E(var|varhat)
  post.mixmean = post.IG.paramb/(post.IG.parama-1)
  post.mean = apply(post.pi*post.mixmean,1,sum)
  
  return(list(pi=post.pi,mean=post.mean,gammaa=post.IG.parama,gammab=post.IG.paramb))
}

# Log-likelihood: L(varhat|c,pi,alpha.vec)
loglike = function(logc,N,K,alpha.vec,n,varhat,pi,unimodal){
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  
  c=exp(logc)
  pimat = exp(outer(rep(1,N),log(pi))+n/2*log(n/2)-lgamma(n/2)
                  +(n/2-1)*outer(log(varhat),rep(1,K))
                  +outer(rep(1,N),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+n/2))
                  -outer(rep(1,N),alpha.vec+n/2)*log(outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  #classprob=pimat/rowSums(pimat)
  logl = sum(log(rowSums(pimat)))
  return(-logl)
}

# Gradient of funtion loglike (w.r.t logc)
gradloglik = function(logc,N,K,alpha.vec,n,varhat,pi,unimodal){
  c=exp(logc)
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  pimat = exp(outer(rep(1,N),log(pi))+n/2*log(n/2)-lgamma(n/2)
            +(n/2-1)*outer(log(varhat),rep(1,K))
            +outer(rep(1,N),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+n/2))
            -outer(rep(1,N),alpha.vec+n/2)*log(outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  classprob = pimat/rowSums(pimat)
  gradmat = c*classprob*(outer(rep(1,N),alpha.vec/c)
                     -outer(rep(1,N),(alpha.vec+n/2)*modalpha.vec)/
    (outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  grad = sum(-gradmat)
  return(grad)
}

# Gradient of funtion loglike for single component prior.se (w.r.t logalpha)
gradloglik.a = function(logc,N,K,logalpha.vec,n,varhat,pi,unimodal){
  alpha.vec=exp(logalpha.vec)
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  c=exp(logc)
  grad=-alpha.vec*sum(logc+log(modalpha.vec)+alpha.vec/modalpha.vec-digamma(alpha.vec)+digamma(alpha.vec+n/2)
           -c*(alpha.vec+n/2)/(c*modalpha.vec+n/2*varhat)-log(c*modalpha.vec+n/2*varhat))
  return(grad)
}


# EM algorithm to estimate pi.se (mixture proportion for varhat's prior),
# c, and alpha(only when singlecomp==TRUE)
# prior.se: nullbiased: add weight to the null component at each iteration; (NEED TESTING!)
#        uniform: uniform weight
# ltol: tolerance of convergence
# maxiter: max number of iterations
EMest_semix = function(varhat,n,c,alpha.vec,pi,prior.se='uniform', ltol=0.0001, maxiter=5000, unimodal, SGD, singlecomp){
  K = length(pi)
  N = length(varhat)
  nullcomp = which.max(alpha.vec)
  if(is.null(c)){
    est.c = TRUE
    c=mean(varhat)
  }else{
    est.c = FALSE
  }
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
    
  if(prior.se=="nullbiased"){ 
    prior.se = rep(1,K)
    prior.se[nullcomp] = 10
  }else if(prior.se=="uniform"){
    prior.se = rep(1,K)
  }else if(prior.se=='unipenalty'){
    prior.se=rep(10,K)
  }
  loglik = rep(NA,maxiter)
  
  mm = post_pi(n,varhat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  loglik[1] = sum(log(m.rowsum))
  classprob = mm/m.rowsum
  
  logc=log(c)
  for(i in 2:maxiter){    
    if (est.c==TRUE & SGD==FALSE){
      est=nlminb(logc,loglike,gradloglik,N=N,K=K,alpha.vec=alpha.vec,n=n,varhat=varhat,pi=pi,unimodal=unimodal)
      c=exp(est$par[1])  
    }else if (est.c==TRUE & SGD==TRUE){
      logc=log(c)-0.001/sqrt(i)*gradloglik(log(c),N,K,alpha.vec,n,varhat,pi,unimodal)
      c=exp(logc)
    } 
    if (singlecomp==TRUE){
      logalpha.vec=max(log(2),log(alpha.vec)-0.001/(i^0.25)*gradloglik.a(log(c),N,K,log(alpha.vec),n,varhat,pi,unimodal))
      alpha.vec=exp(logalpha.vec)
      if(unimodal=='variance'){
        modalpha.vec=alpha.vec+1
      }else if(unimodal=='precision'){
        modalpha.vec=alpha.vec-1
      }
    }
    pi = colSums(classprob)+prior.se-1
    #pi=colSums(classprob)
    pi = ifelse(pi<0,0,pi); pi=pi/sum(pi);
    mm = post_pi(n,varhat,alpha.vec,modalpha.vec,c,pi)
    m.rowsum = rowSums(mm)
    loglik[i] = sum(log(m.rowsum))
    classprob = mm/m.rowsum
    if(abs(loglik[i]-loglik[i-1])<ltol) break;      
  }
  converged = (i< maxiter)
  niter = min(c(i,maxiter))
  return(list(pi=pi,classprob=classprob,loglik.final=loglik,converged=converged,niter=niter,c=c,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec))  
}

logABF = function(betahat,sebetahat,sigmaa){
  T = betahat/sebetahat
  lambda = sebetahat^2/(sebetahat^2+sigmaa^2)
  return(0.5*log(lambda) + 0.5*T^2 *(1-lambda))
  # the following line is same up to a constant, and probably faster:
  # return(dnorm(betahat,0,sqrt(sebetahat^2+sigmaa^2),log=TRUE))
}

#return matrix of ABFs for vector of sigma-a values
#normalized by maximum of each column
#betahat is n vector, sebetahat is n vector, sigmaavec is k vector
#return is n by k matrix of ABFs
matrixABF = function(betahat, sebetahat, sigmaavec){
  k = length(sigmaavec)
  n = length(betahat)
  labf = matrix(0,nrow=n, ncol=k)
  for(i in 1:k){
    labf[,i] = logABF(betahat,sebetahat,sigmaavec[i])
  }
  maxlabf = apply(labf, 1, max)
  labf = labf - maxlabf
  return(exp(labf))
}

totalABF = function(betahat, pi.lik, sigma.lik, sigma.prior){
  L = length(sigma.prior)
  K = dim(sigma.lik)[2]
  n = length(betahat)
  
  abf = 0
  for (j in 1:K){
    labf = matrix(0,nrow=n, ncol=L)
    for(i in 1:L){
      labf[,i] = logABF(betahat,sigma.lik[,j],sigma.prior[i])
    }
    maxlabf = apply(labf, 1, max)
    labf = labf - maxlabf
    abf = abf + outer(pi.lik[,j],rep(1,L))*exp(labf)
  }
  return(abf)
}


# EM algorithm to estimate mixture proportion for beta's prior
# prior.beta: nullbiased: add weight to the null component at each iteration; (NEED TESTING!)
#        uniform: uniform weight
# ltol: tolerance of convergence
# maxiter: max number of iterations
EMest_meanmix = function(betahat, pi.lik, sigma.lik, sigma.prior, prior.beta='nullbiased', ltol=0.0001, maxiter=5000){
  K = length(sigma.lik)
  L = length(sigma.prior)
  N = length(betahat)
  nullcomp = which.min(sigma.prior)
  
  if(prior.beta=="nullbiased"){ 
    prior.beta = rep(1,L)
    prior.beta[nullcomp] = 10
  }else if(prior.se=="uniform"){
    prior.beta = rep(1,L)
  }
  pi = rep(1,L)
  pi[nullcomp] = L
  pi = pi/sum(pi)
  
  loglik = rep(NA,maxiter)
  abf = totalABF(betahat,pi.lik,sigma.lik,sigma.prior)
  m  = t(pi * t(abf)) # abf is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik[1] = sum(log(m.rowsum))
  classprob = m/m.rowsum #an n by k matrix
  
  for(i in 2:maxiter){
    pi = colSums(classprob)+prior.beta-1
    pi = ifelse(pi<0,0,pi); pi=pi/sum(pi);
    abf = totalABF(betahat,pi.lik,sigma.lik,sigma.prior)
    m  = t(pi * t(abf)) # abf is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    loglik[i] = sum(log(m.rowsum))
    classprob = m/m.rowsum
    if(abs(loglik[i]-loglik[i-1])<ltol) break;      
  }
  converged = (i< maxiter)
  niter = min(c(i,maxiter))
  return(list(pi=pi,classprob=classprob,loglik.final=loglik,converged=converged,niter=niter))  
}

#return posterior of being <T (or >T) for a mixture of Gaussians
# each of pi1, mu1, sigma1 is a N by K matrix
# jth column provides parameters for jth mixture of gauusians 
# return an n vector of probabilities
pnormmix = function(T,pi1,mu1,sigma1,lower.tail=TRUE){
  return(apply(pi1 * pnorm(T,mu1,sigma1,lower.tail),1,sum))
}

#return the posterior on beta given a prior
#that is a mixture of normals (pi0,mu0,sigma0)
#and observation betahat \sim N(beta,sebetahat)
#current ABF is only for mu0=0, so would need to
#generalize that for general application
#INPUT: priors: pi0, mu0, sigma0, all k vectors
#       data, betahat (n vector), sebetahat (n vector)
#OUTPUT list (pi1,mu1,sigma1) whose components are each k by n matrices
#k is number of mixture components, n is number of observations
old_posterior_dist = function(pi0,mu0,sigma0,betahat,sebetahat){
  k= length(pi0)
  n= length(betahat)
  
  pi1 = pi0 * t(matrixABF(betahat,sebetahat,sigma0))
  pi1 = apply(pi1, 2, normalize) #pi1 is now an k by n matrix
  
  #make k by n matrix versions of sigma0^2 and sebetahat^2
  # and mu0 and betahat
  s0m2 = matrix(sigma0^2,nrow=k,ncol=n,byrow=FALSE)
  sebm2 = matrix(sebetahat^2,nrow=k,ncol=n, byrow=TRUE)
  mu0m = matrix(mu0,nrow=k,ncol=n,byrow=FALSE)
  bhatm = matrix(betahat,nrow=k,ncol=n,byrow=TRUE)
  
  sigma1 = (s0m2*sebm2/(s0m2 + sebm2))^(0.5)  
  w = sebm2/(s0m2 + sebm2)
  mu1 = w*mu0m + (1-w)*bhatm
  
  #WHERE DATA ARE MISSING, SET POSTERIOR = PRIOR
  ismiss = (is.na(betahat) | is.na(sebetahat)) 
  pi1[,ismiss] = pi0
  mu1[,ismiss] = mu0
  sigma1[,ismiss] = sigma0
  
  return(list(pi=t(pi1),mu=t(mu1),sigma=t(sigma1)))
}

# Compute posterior distribution P(beta|betahat,varhat)
post_distn_mean = function(pi.lik, pi.prior, betahat, sigma.lik, sigma.prior){
  K = dim(sigma.lik)[2]
  L = length(sigma.prior)
  N = length(betahat)
  postpi = NULL
  postmu = NULL
  postsigma = NULL
  for (i in 1:K){
    temppost = old_posterior_dist(pi.prior,0,sigma.prior,betahat,sigma.lik[,i])
    postpi = cbind(postpi, outer(pi.lik[,i],rep(1,L))*temppost$pi)
    postmu = cbind(postmu, temppost$mu)
    postsigma = cbind(postsigma, temppost$sigma)
  }
  return(list(pi=postpi,mu=postmu,sigma=postsigma))
}

#find point estimates of beta from a posterior produced by post_distn_mean
posterior_mean = function(post){
  return(rowSums(post$pi * post$mu))
}

# Approximate t-mixture likelihood P(betahat|varhat,beta)=pi*(beta+b.vec/alpha.vec*t(2*alpha.vec))
# by normal-mixture, with fixed normal params for each component N(0,sigma^2)
# estimate the normal-mixture proportion fi*N(0,sigma^2)
# alpha.vec, beta.vec, pi: P(var|varhat)=pi*InvGamma(alpha.vec,b.vec)
# Solve fi s.t. several quantiles of normal-mixture approx that of t-mixture
library(limSolve)
approxlik = function(params,K){
  alpha.vec = params[1:K]
  b.vec = params[(K+1):(2*K)]
  pi = params[(2*K+1):(3*K)] 
  sigma = params[(3*K+1):length(params)]
  L = length(sigma)
  
  pts = 2^(seq(-2,2))  # NEED TESTING!
  A = pnorm(outer(pts,rep(1,L)),mean=0,sd=outer(rep(1,L),sigma))
  b = rowSums(outer(rep(1,L),pi)*pt(outer(pts,alpha.vec/b.vec),df=outer(rep(1,L),2*alpha.vec)))
  fi = lsei(A = A, B = b, E = rep(1,L), F = 1, G = diag(L), H = rep(0,L))$X
  
  #xgrid = seq(-10,10,by=0.001)
  #densmat = outer(rep(1,length(xgrid)),alpha.vec/b.vec)*dt(outer(xgrid,alpha.vec/b.vec),df=outer(rep(1,length(xgrid)),2*alpha.vec))
  #truedens = apply(outer(rep(1,length(xgrid)),pi)*densmat,1,sum)
  #densmat = dnorm(outer(xgrid,rep(1,L)),mean=0,sd=outer(rep(1,length(xgrid)),sigma))
  #approxdens = apply(outer(rep(1,length(xgrid)),fi)*densmat,1,sum)
  #return(list(fi=fi,sigma=sigma,xgrid=xgrid,truedens=truedens,approxdens=approxdens))
  return(fi)
}

# Adaptive shrinkage for variances
# Unimodal: 'variance': variances~Mix IG with common mode c
#           'precision': precisions~Mix Gamma with common mode c
# SGD: use stochastic gradient descent to est hyperparams
# singlecomp: fit single component prior.se (est both prior.se mean and var params by EB)
vash = function(betahat,sebetahat,n,localfdr=TRUE,sigma.prior=NULL,prior.se='uniform',prior.beta='nullbiased',unimodal='precision',alpha.vec=NULL,c=NULL,pi.se=NULL,SGD=TRUE,singlecomp=FALSE){
  varhat = sebetahat^2
  N = length(betahat)
  if(is.null(alpha.vec)){
    alpha.vec = 2+2^seq(-3,10)
  }
  if(singlecomp==TRUE){
    alpha.vec=max(1/(mean(varhat)^2*sd(varhat)^2),2)
  }
  if(is.null(pi.se)){
    pi.se = rep(1,length(alpha.vec))/length(alpha.vec)
  }
  if(is.null(sigma.prior)){
    sigma.prior = autoselect.mixsd(betahat,sqrt(varhat),10)
  }
  
  # Estimate hyperparams for var prior P(var)~mix IG
  pifit.se = EMest_semix(varhat,n,c,alpha.vec,pi.se,prior.se, ltol=0.0001, maxiter=5000, unimodal,SGD, singlecomp)
  
  # Posterior distribution P(var|varhat)
  post.se = post_distn_se(n,varhat,pifit.se$alpha.vec,pifit.se$modalpha.vec,pifit.se$c,pifit.se$pi)
  
  # Approximate t-mixture likelihood P(betahat|beta,varhat) by normal mixture 
  sigma.lik = outer(rep(1,N),2^seq(-2,2))*apply(post.se$gammaa/post.se$gammab,1,mean)*n/(n-1.9)
  K = dim(sigma.lik)[2]
  post.se.params = cbind(post.se$gammaa, post.se$gammab, post.se$pi,sigma.lik)
  pi.lik = t(apply(post.se.params,1,approxlik, K=length(pifit.se$pi)))
  
  # Estimate hyperparams for beta prior P(beta)~mix normal
  pifit.beta = EMest_meanmix(betahat, pi.lik, sigma.lik, sigma.prior, prior.beta)
  
  # Compute posterior P(beta|betahat,varhat)
  post.beta = post_distn_mean(pi.lik, pifit.beta$pi, betahat, sigma.lik, sigma.prior)
  PosteriorMean = posterior_mean(post.beta)
  
  # Estimate lfdr, lfsr
  PositiveProb = pnormmix(0,post.beta$pi,post.beta$mu,post.beta$sigma,lower.tail=FALSE)
  ZeroProb = rowSums(post.beta$pi[,rep(sigma.prior,K)==0,drop=FALSE])
  NegativeProb =  1- PositiveProb-ZeroProb
  if(localfdr==TRUE){
    localfsr = ifelse(PositiveProb<NegativeProb,PositiveProb+ZeroProb,NegativeProb+ZeroProb)
    localfdr = 2* localfsr
    qvalue = qval.from.lfdr(localfdr)
  }else{
    localfdr = NULL
    localfsr = NULL
    qvalue = NULL
  }
    
  return(list(PosteriorMean=PosteriorMean, pifit.se=pifit.se, pifit.beta=pifit.beta,
              post.se=post.se, post.beta=post.beta, alpha=alpha.vec, c=pifit.se$c,
              pi.beta=pifit.beta$pi, pi.se=pifit.se$pi, pi.lik=pi.lik,
              unimodal=unimodal, sigma.prior=sigma.prior, sigma.lik=sigma.lik,
              localfdr=localfdr,localfsr=localfsr, qvalue=qvalue))
}

# function to plot the Empirical Bayes prior of variances/precisions
# xmax: plot density on (0,xmax)
vashEBprior=function(vashobj,xmax){
  xgrid=seq(0.0001,xmax,by=0.01)
  if(vashobj$unimodal=='variance'){
    EBprior.var.sep=dgamma(outer(1/xgrid,rep(1,length(vashobj$alpha))),
                           shape=outer(rep(1,length(xgrid)),vashobj$alpha),
                           rate=vashobj$c*outer(rep(1,length(xgrid)),vashobj$alpha+1))*outer(1/xgrid^2,rep(1,length(vashobj$alpha)))
    EBprior.var=rowSums(outer(rep(1,length(xgrid)),vashobj$pi)*EBprior.se.var.sep)
    EBprior.prec.sep=dgamma(outer(xgrid,rep(1,length(vashobj$alpha))),
                            shape=outer(rep(1,length(xgrid)),vashobj$alpha),
                            rate=vashobj$c*outer(rep(1,length(xgrid)),vashobj$alpha+1))                    
    EBprior.prec=rowSums(outer(rep(1,length(xgrid)),vashobj$pi)*EBprior.se.prec.sep)
  }else if (vashobj$unimodal=='precision'){
    EBprior.var.sep=dgamma(outer(1/xgrid,rep(1,length(vashobj$alpha))),
                           shape=outer(rep(1,length(xgrid)),vashobj$alpha),
                           rate=vashobj$c*outer(rep(1,length(xgrid)),vashobj$alpha-1))*outer(1/xgrid^2,rep(1,length(vashobj$alpha)))
    EBprior.var=rowSums(outer(rep(1,length(xgrid)),vashobj$pi)*EBprior.se.var.sep)
    EBprior.prec.sep=dgamma(outer(xgrid,rep(1,length(vashobj$alpha))),
                            shape=outer(rep(1,length(xgrid)),vashobj$alpha),
                            rate=vashobj$c*outer(rep(1,length(xgrid)),vashobj$alpha-1))                   
    EBprior.prec=rowSums(outer(rep(1,length(xgrid)),vashobj$pi)*EBprior.se.prec.sep)
  }
  return(list(xgrid=xgrid,EBprior.var=EBprior.se.var,EBprior.prec=EBprior.prec, 
              EBpriorvar.sep=EBprior.var.sep, EBprior.prec.sep=EBprior.prec.sep))
}

