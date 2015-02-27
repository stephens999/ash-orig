#THis file contains functions that used to be used in ash
#but were replaced by the mix.R functionality

#compute normal density for n-vector x
#at each of k values of mu and sigma
#OUTPUT: k by n matrix of normal densities
matdnorm = function (x, mu, sigma, log=FALSE) 
{
  k=length(mu)
  n=length(x)
  d = matrix(rep(x,rep(k,n)),nrow=k)
  return(matrix(dnorm(d, mu, sigma, log),nrow=k))
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

#return posterior of being <T (or >T) for a mixture of Gaussians
# each of pi1, mu1, sigma1 is a k by n matrix
# jth column provides parameters for jth mixture of gauusians 
# return an n vector of probabilities
pnormmix = function(T,pi1,mu1,sigma1,lower.tail=TRUE){
  return(apply(pi1 * pnorm(T,mu1,sigma1,lower.tail),2,sum))
}

#return posterior of being <T (or >T) for a mixture of Gaussians
# each of pi1, mu1, sigma1 is a k-vector
# jth element provides parameters for jth mixture of gauusians 
# return a probabilities
pnormmix.vec = function(T,pi1,mu1,sigma1,lower.tail=TRUE){
  return(sum(pi1 * pnorm(T,mu1,sigma1,lower.tail)))
}

#density function of a convolution of uniform[a,b] with normal(0,sigma)
#note, we don't need a<b
dconvolve.uninorm=function(x,a,b,sigma){
  return((pnorm((x-a)/sigma)-pnorm((x-b)/sigma))/(b-a))
}

#return density for  convolution of uniform[0,b]
#where b is a vector, at vector of betahat and standard errors
#INPUT b: a k-vector of parameters
#betahat, sebetahat: each n vectors of observations and standard errors
#OUTPUT a n by k matrix of the densities

dconvolve.uninorm.matrix = function(betahat, sebetahat,b){
  k = length(b)
  n = length(betahat)
  ld = matrix(0,nrow=n, ncol=k)
  for(i in 1:k){
    ld[,i] = log(dconvolve.uninorm(betahat,0,b[i],sebetahat))
  }
  maxld = apply(ld, 1, max)
  ld = ld - maxld
  return(exp(ld))
}

#density of mixture of uniforms[a,b]
#INPUT: pi = mixture proportions (k vector), [a,b] k vectors
dunimix=function(x,pi, a,b){
  l = ifelse(a<b,a,b) #lower
  u = ifelse(a<b,b,a) # upper
  return(sum((pi/(u-l))[x<u & x>l]))
}

