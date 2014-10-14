#library(SQUAREM)
#library(gaussquad)
source('added.R')

# If x is a n-column vector, turn it into n by 1 matrix
# If x is a matrix, keep it
tomatrix = function(x){
  if(is.vector(x)){
    x = as.matrix(x)
  }
  return(x)
}



#estimate mixture proportions of beta's prior by EM algorithm
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#nullcheck indicates whether to check whether the loglike exceeds the null
#(may not want to use if prior is used)
#VB provides an approach to estimate the approximate posterior distribution
#of mixture proportions of sigmaa by variational Bayes method
#(use Dirichlet prior and approximate Dirichlet posterior)
EMest_mean = function(betahat,sebetahat,pilik,g,prior,null.comp=1,nullcheck=TRUE, maxiter=5000, df){ 
  
  pi.init = g$pi
  k = ncomp(g)
  n = length(betahat)
  l = dim(sebetahat)[2]
  group=rep(1:k,l)
  tol = min(0.1/n,1e-4) # set convergence criteria to be more stringent for larger samples
  
  matrix_lik_raw = t(compdens_conv_mixlik(g,betahat,sebetahat,df,pilik))
  matrix_lik = t(rowsum(t(matrix_lik_raw),group))
  
  EMfit = mixEM(matrix_lik,prior,pi.init,tol, maxiter)
  
  pi = EMfit$pihat     
  loglik = EMfit$B # actually return log lower bound not log-likelihood! 
  converged = EMfit$converged
  niter = EMfit$niter
  loglik.final = EMfit$B[length(EMfit$B)]
  
  null.loglik = sum(log(matrix_lik[,null.comp]))  
  
  if(nullcheck==TRUE){ 
    if(null.loglik > loglik.final){ #check whether exceeded "null" likelihood where everything is null
      pi=rep(0,k)
      pi[null.comp]=1
      m  = t(pi * t(matrix_lik)) 
      m.rowsum = rowSums(m)
      loglik = sum(log(m.rowsum))
    }
  }
  
  g$pi=pi
  
  return(list(loglik=loglik.final,null.loglik=null.loglik,
              matrix_lik=matrix_lik,converged=converged,g=g))
}


#estimate mixture proportions of se's prior by EM algorithm
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
EMest_se = function(betahat,sebetahat,g,prior.se,maxiter=5000, v,unimodal, SGD, singlecomp){ 
  
  pi.init = g$pi
  k = ncomp(g)
  n = length(betahat)
  tol = min(0.1/n,1e-4) # set convergence criteria to be more stringent for larger samples
  
  if(unimodal=='variance'){
    c.init=g$beta[1]/(g$alpha[1]+1)
  }else if(unimodal=='precision'){
    c.init=g$beta[1]/(g$alpha[1]-1)
  }
  
  EMfit = IGmixEM(sebetahat, v, c.init, g$alpha, pi.init, prior.se, unimodal, SGD, singlecomp, tol, maxiter)
  
  loglik = EMfit$B # actually return log lower bound not log-likelihood! 
  converged = EMfit$converged
  niter = EMfit$niter
  loglik.final = EMfit$B[length(EMfit$B)]
  
  g$pi=EMfit$pihat 
  g$c=EMfit$chat 
  if(singlecomp==TRUE){
    g$alpha=EMfit$alphahat
  }
  if(unimodal=='variance'){
    g$beta=g$c*(g$alpha+1)
  }else if(unimodal=='precision'){
    g$beta=g$c*(g$alpha-1)
  }
  
  return(list(loglik=loglik.final,converged=converged,g=g))
}


IGmixEM = function(sebetahat, v, c.init, alpha.vec, pi.init, prior.se, unimodal, SGD, singlecomp, tol, maxiter){
  q = length(pi.init)
  n = length(sebetahat)
  
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  
  if(singlecomp==FALSE){
    params.init=c(log(c.init),pi.init,0)
    res = squarem(par=params.init,fixptfn=fixpoint_se, objfn=penloglik_se, 
                  n=n,k=q,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,v=v,sebetahat=sebetahat,prior.se=prior.se,SGD=SGD,
                  control=list(maxiter=maxiter,tol=tol))
    return(list(chat = exp(res$par[1]), pihat=res$par[2:(length(res$par)-1)], B=-res$value.objfn, 
                niter = res$iter, converged=res$convergence))
  }else{
    params.init=c(log(c.init),log(alpha.vec),0)
    res = squarem(par=params.init,fixptfn=fixpoint_se_estalpha, objfn=penloglik_se_estalpha, 
                  n=n,k=q,v=v,sebetahat=sebetahat,unimodal=unimodal,prior.se=prior.se,SGD=SGD,
                  control=list(maxiter=maxiter,tol=tol))
    return(list(chat = exp(res$par[1]), pihat=1, B=-res$value.objfn, 
                niter = res$iter, converged=res$convergence,
                alphahat=exp(res$par[2])))
  }
  
}


# Estimate the single inv-gamma prior distn params (moments matching)
# Prior: s^2~IG(a,b)
# sebetahat^2~s^2*Gamma(df/2,df/2)
momentm = function(sebetahat,df){
  n = length(sebetahat)
  e = 2*log(sebetahat)-digamma(df/2)+log(df/2)
  ehat = mean(e)
  a = solve_trigamma(mean((e-ehat)^2*n/(n-1)-trigamma(df/2)))
  b = a*exp(ehat+digamma(df/2)-log(df/2))
  return(list(a=a,b=b))
}

# Solve trigamma(y)=x
solve_trigamma = function(x){
  if(x > 1e7){
    y.new = 1/sqrt(x)
  }else if (x < 1e-6){
    y.new = 1/x
  }else{    
    y.old = 0.5+1/x
    delta = trigamma(y.old)*(1-trigamma(y.old)/x)/psigamma(y.old,deriv=2)
    y.new = y.old+delta
    while(-delta/y.new <= 1e-8){
      y.old = y.new
      delta = trigamma(y.old)*(1-trigamma(y.old)/x)/psigamma(y.old,deriv=2)
      y.new = y.old+delta
    }
  }
  return(y.new)
}

# prior.se of se: se|pi,alpha.vec,c ~ pi*IG(alpha_i,c*(alpha_i-1))
# Likelihood: sebetahat^2|se^2 ~ sj*Gamma(v/2,v/2)
# pi, alpha.vec, c: known
# Posterior weight of P(se|sebetahat) (IG mixture distn)
post_pi_vash = function(v,sebetahat,alpha.vec,modalpha.vec,c,pi){
  n = length(sebetahat)
  k = length(alpha.vec)
  post.pi.mat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                       +(v/2-1)*outer(2*log(sebetahat),rep(1,k))
                                       +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                       -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sebetahat^2,c*modalpha.vec,FUN="+")))
  return(pimat=post.pi.mat)
}

fixpoint_se = function(params,n,k,alpha.vec,modalpha.vec,v,sebetahat,prior.se,SGD){
  logc=params[1]
  pi=params[2:(length(params)-1)]
  iter=params[length(params)]
  
  iter=iter+1
  mm = post_pi_vash(v,sebetahat,alpha.vec,modalpha.vec,exp(logc),pi)
  m.rowsum = rowSums(mm)
  classprob = mm/m.rowsum
  newpi = colSums(classprob)+prior.se-1
  newpi = ifelse(newpi<1e-5,1e-5,newpi)
  newpi = newpi/sum(newpi);
  if (SGD==FALSE){
    est=nlminb(logc,loglike.se,gradloglike.se,n=n,k=k,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,v=v,sebetahat=sebetahat,pi=newpi)
    newc=exp(est$par[1])
    newlogc=est$par[1]
  }else{
    newlogc=logc-1/n/sqrt(iter)*gradloglike.se(logc,n,k,alpha.vec,modalpha.vec,v,sebetahat,newpi)  
  }
  params = c(newlogc,newpi,iter)
  return(params)
}

fixpoint_se_estalpha = function(params,n,k,v,sebetahat,unimodal,prior.se,SGD){ 
  logc=params[1]
  alpha.vec=exp(params[2])
  iter=params[3]
  pi=1
  
  iter=iter+1
  
  if (SGD==FALSE){
    if(unimodal=='variance'){
      modalpha=alpha.vec+1
    }else if(unimodal=='precision'){
      modalpha=alpha.vec-1
    }
    estc=nlminb(logc,loglike.se,gradloglike.se,n=n,k=k,alpha.vec=alpha.vec,modalpha.vec=modalpha,v=v,sebetahat=sebetahat,pi=pi)
    newlogc=estc$par[1]
    estalpha=nlminb(log(alpha.vec),loglike.se.a,gradloglike.se.a,c=exp(newlogc),n=n,k=k,v=v,sebetahat=sebetahat,pi=pi,unimodal=unimodal)
    newlogalpha=estalpha$par[1]
  }else{
    newlogalpha=max(log(1),log(alpha.vec)-1/n/sqrt(iter)*gradloglike.se.a(log(alpha.vec),exp(logc),n,k,v,sebetahat,pi,unimodal))
    newalpha=exp(newlogalpha)
    if(unimodal=='variance'){
      newmodalpha=newalpha+1
    }else if(unimodal=='precision'){
      newmodalpha=newalpha-1
    }
    newlogc=logc-1/n/sqrt(iter)*gradloglike.se(logc,n,k,newalpha,newmodalpha,v,sebetahat,pi) 
  }
  
  params = c(newlogc,newlogalpha,iter)
  return(params)
}

penloglik_se = function(params,n,k,alpha.vec,modalpha.vec,v,sebetahat,prior.se,SGD){
  c=exp(params[1])
  pi=params[2:(length(params)-1)]
  iter=params[length(params)]
  priordens = sum((prior.se-1)*log(pi))
  mm = post_pi_vash(v,sebetahat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  loglik = sum(log(m.rowsum))
  return(-(loglik+priordens))
}

penloglik_se_estalpha = function(params,n,k,v,sebetahat,unimodal,prior.se,SGD){
  c=exp(params[1])
  alpha.vec=exp(params[2])
  pi=1
  iter=params[3]
  
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  
  priordens = sum((prior.se-1)*log(pi))
  mm = post_pi_vash(v,sebetahat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  loglik = sum(log(m.rowsum))
  return(-(loglik+priordens))
}

# Log-likelihood: L(sebetahat^2|c,pi,alpha.vec)
loglike.se = function(logc,n,k,alpha.vec,modalpha.vec,v,sebetahat,pi){  
  c=exp(logc)
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sebetahat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sebetahat^2,c*modalpha.vec,FUN="+")))
  #classprob=pimat/rowSums(pimat)
  logl = sum(log(rowSums(pimat)))
  return(-logl)
}

# Gradient of funtion loglike.se (w.r.t logc)
gradloglike.se = function(logc,n,k,alpha.vec,modalpha.vec,v,sebetahat,pi){
  c=exp(logc)
  
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sebetahat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sebetahat^2,c*modalpha.vec,FUN="+")))
  classprob = pimat/rowSums(pimat)
  gradmat = c*classprob*(outer(rep(1,n),alpha.vec/c)
                         -outer(rep(1,n),(alpha.vec+v/2)*modalpha.vec)/
                           (outer(v/2*sebetahat^2,c*modalpha.vec,FUN='+')))
  grad = sum(-gradmat)
  return(grad)
}

# Log-likelihood: L(sebetahat|c,pi,alpha.vec)
loglike.se.a = function(logalpha.vec,c,n,k,v,sebetahat,pi,unimodal){  
  alpha.vec=exp(logalpha.vec)
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sebetahat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sebetahat^2,c*modalpha.vec,FUN="+")))
  #classprob=pimat/rowSums(pimat)
  logl = sum(log(rowSums(pimat)))
  return(-logl)
}

# Gradient of funtion loglike.se for single component prior.se (w.r.t logalpha)
gradloglike.se.a = function(logalpha.vec,c,n,k,v,sebetahat,pi,unimodal){
  alpha.vec=exp(logalpha.vec)
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  grad=-alpha.vec*sum(log(c)+log(modalpha.vec)+alpha.vec/modalpha.vec-digamma(alpha.vec)+digamma(alpha.vec+v/2)
                      -c*(alpha.vec+v/2)/(c*modalpha.vec+v/2*sebetahat^2)-log(c*modalpha.vec+v/2*sebetahat^2))
  return(grad)
}

# Approximate non-standard mixture t-likelihood by normal-mixture
# component i: ~(sqrt(beta[i]/alpha[i]))*T(2*alpha[i]), w.p. pi[i]
approxlik_gq=function(params,q,appsigma,appweight){
  alpha = params[1:q] 
  beta = params[(q+1):(2*q)]
  pi = params[(2*q+1):(3*q)]
  
  fi=numeric(0)
  sigma=numeric(0)
  for (i in 1:q){
    fi = c(fi,pi[i]*appweight[i,])
    sigma = c(sigma, sqrt(beta[i]/alpha[i])*appsigma[i,])
  }
  fi=fi[order(sigma)]
  sigma=sort(sigma)
  return(c(fi,sigma))
}

# Approximate t-distribution (of df) by r-components normal mixture
approxt = function(df, r){
  alpha=df/2-1
  rules=glaguerre.quadrature.rules(r,alpha,normalized=TRUE)
  sigma=sqrt(df/(2*rules[[r]]$x))
  weight=rules[[r]]$w/sum(rules[[r]]$w)
  return(list(sigma=sigma,weight=weight))
}

# Approximate mixture t likelihood (with l components) 
# by mixture normal (with q components)
# pi, alpha, beta are n by l matrices
# component i: ~(sqrt(beta[n,i]/alpha[n,i]))*T(2*alpha[n,i]), w.p. pi[n,i]
mixlik_sd=function(pi,alpha,beta){
  q=dim(pi)[2]
  post.se.params = cbind(alpha, beta, pi)
  ll=max(5,floor(20/q))
  
  appweight=matrix(rep(0,ll*q),nrow=q)
  appsigma=matrix(rep(0,ll*q),nrow=q)
  for (i in 1:q){
    app = approxt(df=2*t(alpha)[i],r=ll) # l components for approximating each t-distribution
    appweight[i,]=app$weight
    appsigma[i,]=app$sigma
  }
  
  results = t(apply(post.se.params,1,approxlik_gq,q=q,appsigma=appsigma,appweight=appweight))
  pilik=results[,1:(dim(results)[2]/2)]
  selik=results[,(dim(results)[2]/2+1):dim(results)[2]]
  return(list(pilik=pilik,selik=selik))
}





#' @useDynLib vash
#todo
#
#' @title Main Adaptive SHrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sebetahat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta.
#'
#' @details See readme for more details
#' 
#' @param betahat, a p vector of estimates 
#' @param sebetahat, a p vector of corresponding standard errors
#' @param method: specifies how ash is to be run. Can be "shrinkage" (if main aim is shrinkage) or "fdr" (if main aim is to assess fdr or fsr)
#' This is simply a convenient way to specify certain combinations of parameters: "shrinkage" sets pointmass=FALSE and prior="uniform";
#' "fdr" sets pointmass=TRUE and prior="nullbiased".
#' @param mixcompdist: distribution of components in mixture ("normal", "uniform" or "halfuniform")
#'
#' @param lambda1: multiplicative "inflation factor" for standard errors (like Genomic Control)
#' @param lambda2: additive "inflation factor" for standard errors (like Genomic Control)
#' @param nullcheck: whether to check that any fitted model exceeds the "null" likelihood
#' in which all weight is on the first component
#' @param df: appropriate degrees of freedom for (t) distribution of betahat/sebetahat
#' @param randomstart: bool, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
#' @param pointmass: bool, indicating whether to use a point mass at zero as one of components for a mixture distribution
#' @param onlylogLR: bool, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, lfdr...
#' @param singlecomp: bool, indicating whether to use a single inverse-gamma distribution as the prior distribution of the variances
#' @param SGD: bool, indicating whether to use the stochastic gradient descent method to fit the prior distribution of the variances
#' @param unimodal: unimodal constraint for the prior distribution of the variances ("variance") or the precisions ("precision")
#' @param prior: string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or 1,1...,1; also can be "nullbiased" 1,1/k-1,...,1/k-1 to put more weight on first component)
#' @param mixsd: vector of sds for underlying mixture components 
#' @param gridmult: the multiplier by which the default grid values for mixsd differ by one another. (Smaller values produce finer grids)
#' @param minimal_output: if TRUE, just outputs the fitted g and the lfsr (useful for very big data sets where memory is an issue) 
#' @param g: the prior distribution for beta (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param maxiter: maximum number of iterations of the EM algorithm
#' 
#'
#' @return a list with elements fitted.g is fitted mixture
#' logLR : logP(D|mle(pi)) - logP(D|null)
#' 
#' @export
#' 
#' @examples 
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#' summary(beta.ash)
#' plot(betahat,beta.ash$PosteriorMean,xlim=c(-4,4),ylim=c(-4,4))
vash = function(betahat,sebetahat,df,method = c("shrink","fdr"), 
                mixcompdist = c("normal","uniform","halfuniform"),
                lambda1=1,lambda2=0,nullcheck=TRUE,randomstart=FALSE, 
                pointmass = FALSE, 
                onlylogLR = FALSE, 
                singlecomp = FALSE,
                SGD = TRUE,
                unimodal = c("variance","precision"),
                prior.mean=c("uniform","nullbiased"), 
                prior.se=NULL,
                mixsd=NULL,gridmult=sqrt(2),
                minimaloutput=FALSE,
                g.se=NULL,
                g.mean=NULL,
                maxiter.mean = 5000,
                maxiter.se = 5000){
  
  #method provides a convenient interface to set a particular combinations of parameters for prior an
  #If method is supplied, use it to set up specific values for these parameters; provide warning if values
  #are also specified by user
  #If method is not supplied use the user-supplied values (or defaults if user does not specify them)
  
  if(!missing(method)){
    method = match.arg(method) 
    if(method=="shrink"){
      if(missing(prior.mean)){
        prior.mean = "uniform"
      } else {
        warning("Specification of prior overrides default for method shrink")
      }
      if(missing(pointmass)){
        pointmass=FALSE
      } else {
        warning("Specification of pointmass overrides default for method shrink")
      }
    }
    
    if(method=="fdr"){
      if(missing(prior.mean)){
        prior.mean = "nullbiased"
      } else {
        warning("Specification of prior overrides default for method fdr")
      }
      if(missing(pointmass)){
        pointmass=TRUE
      } else {
        warning("Specification of pointmass overrides default for method fdr")
      }
    }  
  }
  
  if(missing(unimodal)){
    unimodal = match.arg(unimodal) 
  }
  if(!is.element(unimodal,c("variance","precision"))) stop("Error: invalid type of singlecomp")
  
  if(onlylogLR){
    pointmass = TRUE  
  }
  
  mixcompdist = match.arg(mixcompdist)
  
  if(!is.numeric(prior.mean)){
    prior.mean = match.arg(prior.mean)
  }  
  
  if(length(sebetahat)==1){
    sebetahat = rep(sebetahat,length(betahat))
  }
  if(length(sebetahat) != length(betahat)){
    stop("Error: sebetahat must have length 1, or same length as betahat")
  }
  
  completeobs = (!is.na(betahat) & !is.na(sebetahat))
  n=sum(completeobs)
  
  if(n==0){
    if(onlylogLR){
      return(list(pi=NULL, logLR = 0))
    }
    else{
      stop("Error: all input values are missing")
    }
  }  
  
  if(!missing(g.se)){
    maxiter.se = 1 # if g is specified, don't iterate the EM
    prior.se = rep(1,ncomp(g.se)) #prior is not actually used if g specified, but required to make sure EM doesn't produce warning
    l = ncomp(g.se)
  } else {   
    mm = momentm(sebetahat,df)
    if(singlecomp==TRUE){
      alpha = mm$a
    }else{
      if(mm$a>1){
        alpha = 1+2^seq(-3,3)*(mm$a-1)
      }else{
        alpha = mm$a*2^seq(0,6)
      }
    }
    
    if(unimodal=='precision'){
      alpha = unique(pmax(alpha,1+1e-5)) # alpha<=1 not allowed
      beta = mm$b/(mm$a-1)*(alpha-1)
    }else if(unimodal=='variance'){
      beta = mm$b/(mm$a+1)*(alpha+1)
    }
    
    l = length(alpha)
    
    if(missing(prior.se)){
      prior.se = rep(1,l)
    }
    
    if(randomstart){
      pi.se = rgamma(l,1,1)
    } else {   
      pi.se=rep(1,l)/l
    }
    pi.se=pi.se/sum(pi.se)
    
    g.se=igmix(pi.se,alpha,beta)
  }
  
  if(length(prior.se)!=l){
    stop("invalid prior specification")
  }
  
  pi.fit.se = EMest_se(betahat,sebetahat,g.se,prior.se,maxiter.se,df,unimodal, SGD, singlecomp)
  post.se = post.igmix(pi.fit.se$g,betahat[completeobs],sebetahat[completeobs],df)
  postpi.se = t(comppostprob(pi.fit.se$g,betahat[completeobs],sebetahat[completeobs],df))
  
  PosteriorMean.se = rep(0,length=n)
  PosteriorSD.se = rep(0,length=n)
  PosteriorMean.se[completeobs] = sqrt(postmean(pi.fit.se$g,betahat[completeobs],sebetahat[completeobs],df))
  #PosteriorSD.se[completeobs] = postsd(pi.fit.se$g,betahat[completeobs],sebetahat[completeobs],df) 
  
  if(mixcompdist=="normal"){
    appnorm = mixlik_sd(postpi.se,post.se$alpha,post.se$beta)
    pilik = appnorm$pilik
    selik = appnorm$selik
    moddf = NULL
  }else if(mixcompdist=="uniform" | mixcompdist=="halfuniform"){
    pilik = postpi.se
    selik = sqrt(post.se$beta/post.se$alpha)
    moddf = 2*post.se$alpha[1,]
  }
  pilik = tomatrix(pilik)
  selik = tomatrix(selik)
  l = dim(pilik)[2]
  
  
  if(!missing(g.mean)){
    maxiter.mean = 1 # if g is specified, don't iterate the EM
    prior.mean = rep(1,ncomp(g.mean)) #prior is not actually used if g specified, but required to make sure EM doesn't produce warning
    null.comp=1 #null.comp also not used, but required 
  } else {
    if(is.null(mixsd)){
      mixsd = autoselect.mixsd(betahat[completeobs],sebetahat[completeobs],gridmult)
    }
    if(pointmass){
      mixsd = c(0,mixsd)
    }
    
    null.comp = which.min(mixsd) #which component is the "null"
    
    k = length(mixsd)
    if(!is.numeric(prior.mean)){
      if(prior.mean=="nullbiased"){ # set up prior to favour "null"
        prior.mean = rep(1,k)
        prior.mean[null.comp] = 10 #prior 10-1 in favour of null
      }else if(prior.mean=="uniform"){
        prior.mean = rep(1,k)
      }
    }
    
    if(length(prior.mean)!=k | !is.numeric(prior.mean)){
      stop("invalid prior specification")
    }
    
    if(randomstart){
      pi.mean = rgamma(k,1,1)
    } else {
      if(k<n){
        pi.mean=rep(1,k)/n #default initialization strongly favours null; puts weight 1/n on everything except null
        pi.mean[null.comp] = (n-k+1)/n #the motivation is data can quickly drive away from null, but tend to drive only slowly toward null.
      } else {
        pi.mean=rep(1,k)/k
      }
    }
    
    pi.mean=pi.mean/sum(pi.mean)
    if(!is.element(mixcompdist,c("normal","uniform","halfuniform"))) stop("Error: invalid type of mixcompdist")
    if(mixcompdist=="normal") g.mean=normalmix(pi.mean,rep(0,k),mixsd)
    if(mixcompdist=="uniform") g.mean=unimix(pi.mean,-mixsd,mixsd)
    if(mixcompdist=="halfuniform"){
      g.mean = unimix(c(pi,pi)/2,c(-mixsd,rep(0,k)),c(rep(0,k),mixsd))
      prior.mean = rep(prior.mean, 2)
      pi.mean = rep(pi, 2)
    }
  }
  
  
  pi.fit=EMest_mean(betahat[completeobs],lambda1*selik+lambda2,pilik,g.mean,prior.mean,null.comp=null.comp,nullcheck=nullcheck,maxiter = maxiter.mean, df=moddf)  
  
  if(onlylogLR){
    logLR = tail(pi.fit$loglik,1) - pi.fit$null.loglik
    return(list(fitted.g=pi.fit$g, logLR = logLR))
  } else if(minimaloutput){
    n=length(betahat)
    ZeroProb = rep(0,length=n)
    NegativeProb = rep(0,length=n)
    
    #print("normal likelihood")
    ZeroProb[completeobs] = colSums(comppostprob_mixlik(pi.fit$g,betahat[completeobs],sebetahat[completeobs,],moddf,pilik[completeobs,])[comp_sd(pi.fit$g)==0,,drop=FALSE])     
    NegativeProb[completeobs] = cdf_post_mixlik(pi.fit$g, 0, betahat[completeobs],sebetahat[completeobs,],moddf,pilik[completeobs,]) - ZeroProb[completeobs]
    ZeroProb[!completeobs] = sum(mixprop(pi.fit$g)[comp_sd(pi.fit$g)==0])
    NegativeProb[!completeobs] = mixcdf(pi.fit$g,0) 
    
    lfsr = compute_lfsr(NegativeProb,ZeroProb)
    result = list(fitted.g=pi.fit$g,lfsr=lfsr,fit=pi.fit)
    return(result) 
  } else{
    
    
    #     post = posterior_dist(pi.fit$g,betahat,sebetahat)
    n=length(betahat)
    ZeroProb = rep(0,length=n)
    NegativeProb = rep(0,length=n)
    PosteriorMean = rep(0,length=n)
    PosteriorSD=rep(0,length=n)

    pilikco = tomatrix(pilik[completeobs,])
    selikco = tomatrix(selik[completeobs,])
    
    ZeroProb[completeobs] = colSums(comppostprob_mixlik(pi.fit$g, betahat[completeobs], selikco, moddf,  pilikco)[comp_sd(pi.fit$g)==0,,drop=FALSE])    
    NegativeProb[completeobs] = cdf_post_mixlik(pi.fit$g, 0, betahat[completeobs],selikco,moddf, pilikco) - ZeroProb[completeobs]
    PosteriorMean[completeobs] = postmean_mixlik(pi.fit$g,betahat[completeobs],selikco,moddf, pilikco)
    PosteriorSD[completeobs] =postsd_mixlik(pi.fit$g,betahat[completeobs],selikco,moddf, pilikco) 
    
    #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
    ZeroProb[!completeobs] = sum(mixprop(pi.fit$g)[comp_sd(pi.fit$g)==0])
    NegativeProb[!completeobs] = mixcdf(pi.fit$g,0) 
    PosteriorMean[!completeobs] = mixmean(pi.fit$g)
    PosteriorSD[!completeobs] =mixsd(pi.fit$g)  
    PositiveProb =  1- NegativeProb-ZeroProb    
    
    lfsr = compute_lfsr(NegativeProb,ZeroProb)
    lfsra =  compute_lfsra(PositiveProb,NegativeProb,ZeroProb) 
    
    lfdr = ZeroProb
    qvalue = qval.from.lfdr(lfdr)
    
    result = list(fitted.g=pi.fit$g,fitted.g.se=pi.fit.se$g,logLR =tail(pi.fit$loglik,1) - pi.fit$null.loglik,
                  PosteriorMean = PosteriorMean,PosteriorSD=PosteriorSD,PosteriorMean.se=PosteriorMean.se,
                  #PosteriorSD.se=PosteriorSD.se,
                  PositiveProb =PositiveProb,NegativeProb=NegativeProb, ZeroProb=ZeroProb,lfsr = lfsr,lfsra=lfsra, lfdr=lfdr,qvalue=qvalue,
                  fit=pi.fit,fit.se=pi.fit.se,lambda1=lambda1,lambda2=lambda2,call=match.call(),data=list(betahat = betahat, sebetahat=sebetahat))
    class(result)= "ash"
    return(result)
    
  }
  
}

