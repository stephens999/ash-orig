#return ABF for vector of betahat and standard errors
ABF = function(betahat, sebetahat,sigmaa){
  T = betahat/sebetahat
  lambda = sebetahat^2/(sebetahat^2+sigmaa^2)
  return((sqrt(lambda) * exp(0.5*T^2 *(1-lambda))))
}

logABF = function(betahat,sebetahat,sigmaa){
  T = betahat/sebetahat
  lambda = sebetahat^2/(sebetahat^2+sigmaa^2)
  return(0.5*log(lambda) + 0.5*T^2 *(1-lambda))
# the following line is same up to a constant, and probably faster:
# return(dnorm(betahat,0,sqrt(sebetahat^2+sigmaa^2),log=TRUE))
}

#compute normal density for n-vector x
#at each of k values of mu and sigma
#OUTPUT: k by n matrix of normal densities
matdnorm = function (x, mu, sigma) 
{
  k=length(mu)
  n=length(x)
  d = matrix(rep(x,rep(k,n)),nrow=k)
  return(matrix(dnorm(d, mu, sigma),nrow=k))
}

#compute density of mixture of k normals
#for n-vector x
#INPUT: pi, mu and sigma are k vectors; x is n-vector
#OUTPUT: n-vector of densities
mixdnorm = function (x, pi, mu, sigma) 
{
  return (pi %*% matdnorm(x,mu,sigma))
}


#compute mixture log likelihood for data x, for mixture of k normals
#INPUT: x an n vector of data
#mixture component parameters pi, mu, sigma, each being k vectors
#
mixLoglik = function(x,pi,mu,sigma){
  dd=matdnorm(x,mu,sigma) 
  return(sum(log(pi %*% dd)))
}

#compute mixture log likelihood for data x, for mixture of k normals
#plus observation-specific errors with observation-specific variances
#ie x_s is sum_k pi_k N(x_s; mu_k, sigma^2_k + sigma^2_s)
# [optionally, use sigma^2_k * sigma^2_s if FUN="*"]
#INPUT: x an n vector of data
#mixture component parameters pi, mu, sigma, each being k vectors
#se, the observation-specific standard errrors
mixseLoglik = function(x,pi,mu,sigma,se,FUN="+"){
  k=length(mu)
  n=length(x)
  d = matrix(rep(x,rep(k,n)),nrow=k)
  s = sqrt(outer(sigma^2,se^2, FUN)) # k by n matrix of standard errors
  dd = matrix(dnorm(d, mu, s),nrow=k) 
  return(sum(log(pi %*% dd)))
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

#return the KL-divergence between 2 dirichlet distributions
#p,q are the vectors of dirichlet parameters of same lengths
diriKL = function(p,q){
  p.sum = sum(p)
  q.sum = sum(q)
  k = length(q)
  KL = lgamma(q.sum)-lgamma(p.sum)+sum((q-p)*(digamma(q)-digamma(rep(q.sum,k))))+sum(lgamma(p)-lgamma(q))
  return(KL)
}

#estimate mixture proportions of sigmaa by EM algorithm
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#nullcheck indicates whether to check whether the loglike exceeds the null
#(may not want to use if prior is used)
#Introduced sigma.est with intention of implenting an option to esitmate
#sigma rather than fixing it, but this not yet implemented.
#VB provides an approach to estimate the approximate posterior distribution
#of mixture proportions of sigmaa by variational Bayes method
#(use Dirichlet prior and approximate Dirichlet posterior)

EMest = function(betahat,sebetahat,sigmaavec,pi,sigma.est=FALSE,nullcheck=TRUE,prior=NULL,nc=NULL,VB=FALSE,ltol=0.0001, maxiter=5000){ 
  if(!is.null(nc)&sigma.est==TRUE){
    sigmaavec=2^(seq(-15,3,length.out=nc))
  }
  k=length(sigmaavec)
  n = length(betahat)
  sigmamin=min(sigmaavec)
  sigmamax=max(sigmaavec)
  null.comp = which.min(sigmaavec) #which component is the "null"
  if(is.null(prior)){ # set up prior to be 1,1/(k-1),...,1/(k-1) to favour "null"
    prior = rep(1/(k-1),k)
    prior[null.comp] = 1
  }else if(prior=="uniform"){
    prior = rep(1,k)
  }

  abf = matrixABF(betahat,sebetahat,sigmaavec)

  if(VB==TRUE){
    F1 = rep(0,maxiter)
    pipost = prior # Dirichlet posterior on pi
    dnorm.mat = matrix(dnorm(rep(betahat,each=k), rep(0,k*n), sd=sqrt(rep(sebetahat^2,each=k)+rep(sigmaavec^2,n))),ncol=k,byrow=TRUE)
    avgpipost = matrix(exp(digamma(rep(pipost,n))-digamma(rep(sum(pipost),k*n))),ncol=k,byrow=TRUE)
    classprob = avgpipost * dnorm.mat
    classprob = classprob/rowSums(classprob) # n by k matrix
    
    priorSum = sum(prior)
    pipostSum = sum(pipost)
    F1[1] = sum(classprob*log(avgpipost*dnorm.mat)) - diriKL(prior,pipost) #negative free energy
    
    for(i in 2:maxiter){  
      pipost = colSums(classprob) + prior
      
      #Now re-estimate pipost
      avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
      classprob = avgpipost*dnorm.mat
      classprob = classprob/rowSums(classprob) # n by k matrix
      
      priorSum = sum(prior)
      pipostSum = sum(pipost)
      F1[i] = sum(classprob*log(avgpipost*dnorm.mat)) - diriKL(prior,pipost)
      
      if(abs(F1[i]-F1[i-1])<ltol) break;
    }  
    pi = pipost/sum(pipost) # use posterior mean to estimate pi
    
    m  = t(pi * t(abf)) # abf is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    loglik.final = sum(log(m.rowsum))
    null.loglik = sum(log(abf[,null.comp]))
    loglik = F1 # actually return F1 not log-likelihood! 
  }else{
    loglik = rep(0,maxiter)
    m  = t(pi * t(abf)) # abf is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    loglik[1] = sum(log(m.rowsum))
    classprob = m/m.rowsum #an n by k matrix
    
    for(i in 2:maxiter){  
      pi = colSums(classprob) + prior-1
      pi = ifelse(pi<0,0,pi) #set any estimates that are less than zero, which can happen with prior<1, to 0
      pi = normalize(pi)
      
      #estimate sigma
      if(sigma.est==TRUE){
        for(j in 1:k){
          pj=classprob[,j]
          f=function(x) sum((betahat^2/(sebetahat^2+x)^2-1/(sebetahat^2+x))*pj)
          if(f(sigmamin^2)<=0){
            sigmaavec[j]=sigmamin
          }else if(f(sigmamax^2)>=0){
            sigmaavec[j]=sigmamax
          }else{
            sigmaavec[j]=sqrt(uniroot(f,c(sigmamin^2,sigmamax^2))$root)          
          }
        }
        abf = matrixABF(betahat,sebetahat,sigmaavec)
      }
      
      #Now re-estimate pi
      m  = t(pi * t(abf)) 
      m.rowsum = rowSums(m)
      loglik[i] = sum(log(m.rowsum))
      classprob = m/m.rowsum
      
      if(abs(loglik[i]-loglik[i-1])<ltol) break;
    }
    null.loglik = sum(log(abf[,null.comp]))
    loglik.final = loglik[i]
  }
  
  if(nullcheck==TRUE){
    if(null.loglik > loglik.final){ #check whether exceeded "null" likelihood where everything is null
      pi=rep(0,k)
      pi[null.comp]=1
      m  = t(pi * t(abf)) 
      m.rowsum = rowSums(m)
      loglik[i] = sum(log(m.rowsum))
      classprob = m/m.rowsum
    }
  }

  return(list(pi=pi,classprob=classprob,sigmaavec=sigmaavec,loglik=loglik[1:i],null.loglik=null.loglik,
            abf=abf,niter=i, converged = (i<maxiter), temp1=sum(log(abf[,null.comp])),temp2=loglik.final))
}

normalize = function(x){return(x/sum(x))}

#return the posterior on beta given a prior
#that is a mixture of normals (pi0,mu0,sigma0)
#and observation betahat \sim N(beta,sebetahat)
#current ABF is only for mu0=0, so would need to
#generalize that for general application
#INPUT: priors: pi0, mu0, sigma0, all k vectors
#       data, betahat (n vector), sebetahat (n vector)
#OUTPUT list (pi1,mu1,sigma1) whose components are each k by n matrices
#k is number of mixture components, n is number of observations
posterior_dist = function(pi0,mu0,sigma0,betahat,sebetahat){
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
  
  return(list(pi=pi1,mu=mu1,sigma=sigma1))
}



#find mean and variance of a mixture of normals
#INPUT: x is a list with elements pi mu and sigma, each k by n matrices
#OUTPUT; the n vectors of mean and variances of mixtures correspondign to columns of pi, mu and sigma
normmix.mv=function(x){
  Ex = colSums(x$pi * x$mu)
  Ex2 = colSums(x$pi* (x$mu^2 + x$sigma^2))
  Vx = Ex2- Ex^2
  return(list(mean = Ex, var=Vx))
}  

#helper function for posterior_sample
#samples nsamp integers from 1:k according to a given prob vector
sample_component=function(p,nsamp){
  
  return(sample(length(p),nsamp,replace=TRUE,prob=p))
}

#m is a k by n matrix
#comp is a n vector of values in 1-k
#returns the comp[i]-th row of m[,i] 
extract_component=function(comp,m){
  return(m[cbind(comp,seq(comp))])
}

#returns matrix of nsamp samples from posterior
#computed using posterior_dist
# NOTE THIS IS UNTESTED, AND PROBABLY NOT WORKING YET...
posterior_sample = function(post,nsamp){
  component = as.vector(apply(post$pi,2,sample_component,nsamp=nsamp))
  k = ncol(post$pi)
  s = rep(1:k,rep(nsamp,k))
  index = cbind(component,s) #set up indices of mu and sigma to extract
  m = post$mu[index]
  ss = post$sigma[index]
  res = matrix(rnorm(length(m),mean=m,sd=ss),nrow=nsamp)
  return(res)
}


#find point estimates of beta from a posterior produced by posterior_dist
posterior_mean = function(post){
  return(colSums(post$pi * post$mu))
}

#return posterior of being >T for a mixture of Gaussians
# each of pi1, mu1, sigma1 is a k by n matrix
# jth column provides parameters for jth mixture of gauusians 
# return an n vector of probabilities
PosteriorProbExceedsT = function(pi1,mu1,sigma1,T=0){
  return(apply(pi1 * pnorm(T,mu1,sigma1,lower.tail=FALSE),2,sum))
}
  
#return the "effective" estimate
#that is the effect size betanew whose z score betanew/se
#would give the same p value as betahat/se compared to a t with df
effective.effect=function(betahat,se,df){
  p = pt(betahat/se,df)
  qnorm(p,sd=se)
}

get_loglik = function(z.ash){
  return(rev(z.ash$fit$loglik)[1])
}

#compute corresponding q values from a vector of local fdr estimates
#INPUT: localfdr a vector of local fdr estimates
#OUTPUT: qvalue, a vector of q value estimates
qval.from.localfdr = function(localfdr){
  o = order(localfdr)
  qvalue=rep(NA,length(localfdr))
  qvalue[o] = (cumsum(sort(localfdr))/(1:sum(!is.na(localfdr))))
  return(qvalue)
}

# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
autoselect.sigmaavec = function(betahat,sebetahat){
  sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  npoint = ceiling(log2(sigmaamax/sigmaamin))
  return(2^((-npoint):0) * sigmaamax)
}



#main adaptive shrinkage function
#takes a vector of betahats and ses;
#fits a mixture of normals to it
# and returns posteriors
#INPUT: betahat (p vector); sebetahat (p vector of standard errors)
#df: degrees of freedome used to compute sebetahat
#randomstart: bool, indicating whether to initialize EM randomly
#usePointMass: bool, indicating whether to use a point mass at zero as one of components for a mixture distribution
#onlylogLR (= FALSE) : bool, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, localfdr...
#localfdr (=TRUE) : bool,  indicating whether to compute localfdr and q-value
#auto (=FALSE): bool, whether to try to select the sigmaavec vector automatically (beta functionality)
#OUTPUT: 
#logLR : logP(D|mle(pi)) - logP(D|null)
#Things to do: automate choice of sigmavec
# check sampling routin
# check number of iterations
ash = function(betahat,sebetahat,nullcheck=TRUE,df=NULL,randomstart=FALSE, usePointMass = FALSE, onlylogLR = FALSE, localfdr = TRUE, prior=NULL, sigmaavec=NULL, auto=FALSE, sigma.est=FALSE, nc=NULL, VB=FALSE){
  #if df specified, convert betahat so that bethata/sebetahat gives the same p value
  #from a z test as the original effects would give under a t test with df=df
  if(!is.null(df)){
    betahat = effective.effect(betahat,sebetahat,df)
  }	
  
  if(length(sebetahat)==1){
    sebetahat = rep(sebetahat,length(betahat))
  }
  if(is.null(sigmaavec)){
    sigmaavec = c(0.00025,0.0005,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1.024,2.048,4.096,8.192) 
  }
  
  completeobs = !is.na(betahat) & !is.na(sebetahat)
  if(auto==TRUE){
    sigmaavec= autoselect.sigmaavec(betahat[completeobs],sebetahat[completeobs])
  }
  if(usePointMass){
        sigmaavec = c(0,sigmaavec)
  }
  if(sigma.est==TRUE&!is.null(nc)){
    k=nc
  }else{
    k=length(sigmaavec)
  }
  pi = rep(1, k)
  pi[1]=k
  pi=normalize(pi)
  if(randomstart){pi=rgamma(k,1,1)}
  
  pi.fit=EMest(betahat[completeobs],sebetahat[completeobs],sigmaavec=sigmaavec,pi=pi,sigma.est=sigma.est,prior=prior,nullcheck=nullcheck,nc=nc,VB=VB)  
  if(sigma.est==TRUE){
    sigmaavec=pi.fit$sigmaavec
  }
  if(onlylogLR){
	logLR = pi.fit$temp2 - pi.fit$temp1
	return(list(pi=pi.fit$pi, logLR = logLR))
  }else{
   	post = posterior_dist(pi.fit$pi,0,sigmaavec,betahat,sebetahat)
  	PositiveProb = PosteriorProbExceedsT(post$pi,post$mu,post$sigma,0)
  	PosteriorMean = posterior_mean(post)
  	if(localfdr){
   		localfdr = ifelse(PositiveProb<0.5,PositiveProb,1-PositiveProb)
   		qvalue = qval.from.localfdr(localfdr)
  	}else{
   		localfdr=NULL
   		qvalue=NULL
  	}
  	fitted.f= list(pi=pi.fit$pi,sigma=sigmaavec,mu=rep(0,k))
    result = list(post=post,fitted.f=fitted.f,PosteriorMean = PosteriorMean,PositiveProb =PositiveProb,localfdr=localfdr,qvalue=qvalue,fit=pi.fit)
	  class(result)= "ash"
    return(result)

  }
  #if(nsamp>0){
  #  sample = posterior_sample(post,nsamp)
  #}
}

print.ash =function(a){
  a$fitted.f
}

plot.ash = function(a,xmin,xmax){
  x = seq(xmin,xmax,length=1000)
  y = density(a,x)
  plot(x,y,type="l")
}

#compute the predictive density of an observation
#given the fitted ash object a and the vector se of standard errors
#not implemented yet
predictive=function(a,se){
  
}
  
#return the density of the fitted underlying hierarchical f
#a is an ash object
#x is a vector at which the density should be computed
density.ash=function(a,x){mixdnorm(x,a$fitted.f$pi,a$fitted.f$mu,a$fitted.f$sigma)}
