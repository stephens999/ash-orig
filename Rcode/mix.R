############## HELPER FUNCTIONS #################################################

#OUTPUT: k by n matrix of normal densities
matdnorm = function (x, mu, sigma, log=FALSE) 
{
  k=length(mu)
  n=length(x)
  d = matrix(rep(x,rep(k,n)),nrow=k)
  return(matrix(dnorm(d, mu, sigma, log),nrow=k))
}



################################## GENERIC FUNCTIONS ############################
# find matrix of densities at y, for each component of the mixture
# INPUT y is an n-vector
# OUTPUT k by n matrix of densities
matdens = function(x,y,log=FALSE){
  UseMethod("matdens")
}
matdens.default = function(x,y,log=FALSE){
  stop("No such class")
}

#find cdf at y, a generic function
mixcdf = function(x,y,lower.tail=TRUE){
  UseMethod("mixcdf")
}
mixcdf.default = function(x,y,lower.tail=TRUE){
  stop("No such class")
}

#find density at y, a generic function
dens = function(x,y){
  UseMethod("dens")
}
dens.default = function(x,y){
  return (x$pi %*% matdens(x, y))
}

#find log likelihood of data in y
LogLik = function(x,y){
  UseMethod("LogLik")
}
LogLik.default = function(x,y){
  sum(dens(x,y))
}

#compute the density of the components of the mixture m
#when convoluted with a normal with standard deviation s
#the density is evaluated at x
#x and s are n-vectors
#m is a mixture with k components
#output is a k by n matrix of densities
compdens_conv = function(m, x, s){
  UseMethod("compdens_conv")
}
compdens_conv.default = function(m,x, s){
  stop("No such class")
}

#compute density of mixture m convoluted with normal of sd (s)
#at locations x
#m is a mixture
#x is an n vector
#s is an n vector or integer
dens_conv = function(m,x,s){
  UseMethod("dens_conv")
}
dens_conv.default = function(m,x,s){
  colSums(m$pi * compdens_conv(m,x,s))
}

#compute the posterior prob that each observation
#came from each component of the mixture m
#output a k by n vector of probabilities
#computed by weighting the component densities by pi
#and then normalizing
comppostprob=function(m,x,s){
 UseMethod("comppostprob") 
}
comppostprob.default = function(m,x,s){
  t(t(m$pi * compdens_conv(m,x,s))/dens_conv(m,x,s))
}
# evaluate cdf of posterior distribution of beta at c
# m is the prior on beta, a mixture
# c is location of evaluation
# assumption is betahat | beta \sim N(beta,sebetahat)
# m is a mixture with k components
# c a scalar
# betahat, sebetahat are n vectors 
# output is a k by n matrix
compcdf_post=function(m,c,betahat,sebetahat){
  UseMethod("compcdf_post")
}
compcdf_post.default=function(m,c,betahat,sebetahat){
  stop("method compcdf_post not written for this class")
}


cdf_post = function(m,c,betahat,sebetahat){
  UseMethod("cdf_post")
}
cdf_post.default=function(m,c,betahat,sebetahat){
  colSums(comppostprob(m,betahat,sebetahat)*compcdf_post(m,c,betahat,sebetahat))
}


############################### METHODS FOR normalmix class ###########################

# constructor
normalmix = function(pi,mean,sd){
  structure(data.frame(pi,mean,sd),class="normalmix")
}


matdens.normalmix = function(x,y,log=FALSE){
  k=length(x$mean)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dnorm(d, x$mean, x$sd, log),nrow=k))  
}

mixcdf.normalmix = function(x,y,lower.tail=TRUE){
  x$pi %*%  vapply(y,pnorm,x$mean,x$mean,x$sd,lower.tail)
}


############################### METHODS FOR unimix class ###########################

#constructor
#here the uniform is parameterized in terms of min=mean-sd and min=mean+sd
#(so obviously sd is a misnomer!)
unimix = function(pi,a,b){
  structure(data.frame(pi,a,b),class="unimix")
}

matdens.unimix = function(x,y,log=FALSE){
  k=length(x$a)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dunif(d, x$a, x$b, log),nrow=k))  
}

#density of convolution of each component of a unif mixture with N(0,s) at x
# x an n-vector
#return a k by n matrix
compdens_conv.unimix = function(m, x, s){
  return(t(pnorm(outer(x,m$a,FUN="-")/s)
          -pnorm(outer(x,m$b,FUN="-")/s))/(m$b-m$a))
}


#c is a scalar
#m a mixture with k components
#betahat a vector of n observations
#sebetahat an n vector of standard errors
#return a k by n matrix of the posterior cdf
compcdf_post.unimix=function(m,c,betahat,sebetahat){
  c = ifelse(c>m$b,m$b,c)
  c = ifelse(c<m$a,m$a,c)
  pna = pnorm(outer(betahat,m$a,FUN="-")/sebetahat)
  pnc = pnorm(outer(betahat,c,FUN="-")/sebetahat)
  pnb = pnorm(outer(betahat,m$b,FUN="-")/sebetahat)
  t((pnc-pna)/(pnb-pna))
}


################################# TESTING #####################################
temp = normalmix(c(0.5,0.5),c(-3,2),c(1,1))
x = seq(-5,5,length=100)
plot(x,dens(temp,x),type="l")
all.equal(matdens(temp,x),matdnorm(x, temp$mean,temp$sd))

plot(x,mixcdf(temp,x),type="l")

temp2 = unimix(c(0.5,0.5),c(-3,3),c(-1,4))
plot(x,dens(temp2,x),type="l")

plot(x,compdens_conv(temp2,x,0.01)[1,],type="l")
plot(x,compdens_conv(temp2,x,0.1)[2,],type="l")
plot(x,dens_conv(temp2,x,0.2),type="l")

#Note the second of these gives errors - need
#to rewrite to deal with the numerical problems
compcdf_post(temp2,0,c(1,2,3),c(10,10,10))
compcdf_post(temp2,-2,c(-2),c(0.1))
