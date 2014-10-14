#find log likelihood of data, when 
#the mixture m is convolved with l-comp normal mixture
#in betahat with mixture sd betahatsd, mixture proportion pilik
#betahatsd is an n by l matrix
#betahat is an n vector
#v is the degree of freedom
#pilik is an l vector
#' @title loglik_conv_mixlik
#' 
#' @export
#' 
loglik_conv_mixlik = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  UseMethod("loglik_conv_mixlik")
}
#' @title loglik_conv.default
#' 
#' @export
#' 
loglik_conv_mixlik.default = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  sum(log(dens_conv_mixlik(m,betahat,betahatsd,v,pilik,FUN)))
}

#compute the density of the components of the mixture m
#when convoluted with l-components normal mixture with standard deviation s
#or C (scale vector) multiplies scaled (se) student.t l-components mixture with df v
#with mixture proportion pilik
#the density is evaluated at x
#x is an n-vector
#s and pilik are n by l matrices
#v and c are l-vectors
#m is a mixture with k components
#output is a (k*l) by n matrix of densities
compdens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("compdens_conv_mixlik")
}
compdens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  dens=NULL
  for (i in 1:dim(pilik)[2]){
    dens=rbind(dens,pilik[,i]*compdens_conv(m,x,s[,i],v[i],FUN))
  }
  return(dens)
}

#compute density of mixture m convoluted with l-components
#normal mixture of sd (s) or student t mixture with df v
#with mixture proportion pilik
#at locations x
#m is a mixture
#x is an n vector
#s and pilik are n by l matrices
dens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("dens_conv_mixlik")
}
dens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  l=dim(pilik)[2]
  colSums(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik,FUN))
}

#compute the posterior prob that each observation
#came from each component of the mixture m
#output a k by n vector of probabilities
#computed by weighting the component densities by pi
#and then normalizing
#when likelihood is an l-components normal mixture or student t mixture
#with mixture proportion pilik
comppostprob_mixlik=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik") 
}
comppostprob_mixlik.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  group=rep(1:k,l)
  return(rowsum(t(tmp),group))
}

#compute the posterior prob that each observation
#came from each component of the mixture m and the likelihood mixture
#output a (k*l) by n vector of probabilities
#computed by weighting the component densities by pi
#and then normalizing
#when likelihood is an l-components normal mixture or student t mixture
#with mixture proportion pilik
comppostprob_mixlik2=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik2") 
}
comppostprob_mixlik2.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  return(t(tmp))
}

# evaluate cdf of posterior distribution of beta at c
# m is the prior on beta, a mixture
# c is location of evaluation
# assumption is betahat | beta \sim 
# l-components normal mixture with sd sebetahat
# and mixture proportion pilik
# m is a mixture with k components
# c a scalar
# betahat is an n vector 
# sebetahat and pilik are n by l matrices 
# output is a (k*l) by n matrix
compcdf_post_mixlik=function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("compcdf_post_mixlik")
}
compcdf_post_mixlik.default=function(m,c,betahat,sebetahat,v,pilik){
  cdf=NULL
  for (i in 1:dim(pilik)[2]){
    cdf=rbind(cdf,pilik[,i]*compcdf_post(m,c,betahat,sebetahat[,i],v[i]))
  }
  cdf
}

cdf_post_mixlik = function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("cdf_post_mixlik")
}
cdf_post_mixlik.default=function(m,c,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik)*
            compcdf_post_mixlik(m,c,betahat,sebetahat,v,pilik))
}

#output posterior mean for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
postmean_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean_mixlik")
}
postmean_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean_mixlik(m,betahat,sebetahat,v,pilik))
}

#output posterior mean-squared value for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
postmean2_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean2_mixlik")
}
postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean2_mixlik(m,betahat,sebetahat,v,pilik))
}

#output posterior sd for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
postsd_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("postsd_mixlik")
}
postsd_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  sqrt(postmean2_mixlik(m,betahat,sebetahat,v,pilik)-postmean_mixlik(m,betahat,sebetahat,v,pilik)^2)
}

#output posterior mean-squared value for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
comp_postmean2_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean2_mixlik")
}
comp_postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  comp_postsd_mixlik(m,betahat,sebetahat,v,pilik)^2 + 
    comp_postmean_mixlik(m,betahat,sebetahat,v,pilik)^2
}

#output posterior mean for beta for each component of prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
comp_postmean_mixlik=function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean_mixlik")
}
comp_postmean_mixlik.default=function(m,betahat,sebetahat,v,pilik){
  mean=NULL
  for (i in 1:dim(pilik)[2]){
    mean=rbind(mean,comp_postmean(m,betahat,sebetahat[,i],v[i]))
  }
  return(mean)
}

#output posterior sd for beta for each component of prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
comp_postsd_mixlik=function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postsd_mixlik")
}
comp_postsd_mixlik.default=function(m,betahat,sebetahat,v,pilik){
  sd=NULL
  for (i in 1:dim(pilik)[2]){
    sd=rbind(sd,comp_postsd(m,betahat,sebetahat[,i],v[i]))
  }
  return(sd)
}




############################### METHODS FOR igmix class ###########################

#' @title Constructor for igmix class
#'
#' @description Creates an object of class igmix (finite mixture of univariate inverse-gammas)
#' 
#' @details None
#' 
#' @param pi: vector of mixture proportions
#' @param alpha: vector of shape parameters
#' @param beta: vector of rate parameters
#' 
#' @return an object of class igmix
#' 
#' @export
#' 
#' @examples igmix(c(0.5,0.5),c(1,1),c(1,2))
#' 
igmix = function(pi,alpha,beta){
  structure(data.frame(pi,alpha,beta),class="igmix")
}

comp_sd.igmix = function(m){
  m$beta/(m$alpha-1)/sqrt(m$alpha-2)
}

comp_mean.igmix = function(m){
  m$beta/(m$alpha-1)
}

compdens.igmix = function(x,y,log=FALSE){
  k=ncomp(x)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dgamma(1/d, shape=x$alpha, rate=outer(m$beta,1/y^2),log),nrow=k))  
}

#density of product of each component of a inverse-gamma mixture with Gamma(v/2,v/2) at x
# x an n-vector at which density is to be evaluated
#return a k by n matrix
compdens_conv.igmix = function(m,x,s,v,FUN="+"){  
  k=ncomp(m)
  n=length(x)
  dens = t(exp(v/2*log(v/2)-lgamma(v/2)
               +(v/2-1)*outer(log(x^2),rep(1,k))
               +outer(rep(1,n),m$alpha*log(m$beta)-lgamma(m$alpha)+lgamma(m$alpha+v/2))
               -outer(rep(1,n),m$alpha+v/2)*log(outer(v/2*x^2,m$beta,FUN="+"))))
  return(dens)
  
}

#compute posterior shape (alpha1) and rate (beta1)
post.igmix = function(m,betahat,sebetahat,v){
  n = length(sebetahat)
  alpha1 = outer(rep(1,n),m$alpha+v/2)
  beta1 = outer(m$beta,v/2*sebetahat^2,FUN="+")
  ismissing = is.na(sebetahat)
  beta1[,ismissing]=m$beta
  return(list(alpha=alpha1,beta=t(beta1)))
}


comp_cdf.igmix = function(x,y,lower.tail=TRUE){
  vapply(y,pigamma,x$alpha,x$alpha,x$beta,lower.tail)
}

#c is a scalar
#m a mixture with k components
#betahat a vector of n observations
#sebetahat an n vector of standard errors
#return a k by n matrix of the posterior cdf
compcdf_post.igmix=function(m,c,betahat,sebetahat,v){
  #compute posterior shape (alpha1) and rate (beta1)
  alpha1 = m$alpha+v/2
  beta1 = outer(m$beta,v/2*sebetahat^2,FUN="+")
  ismissing = is.na(sebetahat)
  beta1[,ismissing]=m$beta
  
  return(t(pigamma(c,alpha1,beta1)))
}

#return posterior mean for each component of prior m, given observations sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
comp_postmean.igmix = function(m,betahat,sebetahat,v){
  k = length(m$pi)
  n=length(betahat)
  tmp=outer(v/2*sebetahat^2,m$beta,FUN="+")/outer(rep(1,n),m$alpha+v/2-1)
  ismissing = is.na(sebetahat)
  tmp[ismissing,]=m$beta/(m$alpha-1) #return prior mean when missing data
  t(tmp)
}


#return posterior standard deviation for each component of prior m, given observations betahat and sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
comp_postsd.igmix = function(m,betahat,sebetahat,v){
  k = length(m$pi)
  n=length(betahat)
  alpha1 = outer(rep(1,n),m$alpha+v/2-1)
  beta1 = outer(v/2*sebetahat^2,m$beta,FUN="+")
  return(t(beta1/(alpha1-1)/sqrt(alpha1-2)))
}
