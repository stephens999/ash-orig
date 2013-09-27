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



############################### METHODS FOR unimix class ###########################

#constructor
#here the uniform is parameterized in terms of min=mean-sd and min=mean+sd
#(so obviously sd is a misnomer!)
unimix = function(pi,mean,sd){
  structure(data.frame(pi,mean,sd),class="unimix")
}

matdens.unimix = function(x,y,log=FALSE){
  k=length(x$mean)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dunif(d, x$mean-x$sd, x$mean+x$sd, log),nrow=k))  
}



################################# TESTING #####################################
temp = normalmix(c(0.5,0.5),c(-3,2),c(1,1))
x = seq(-5,5,length=100)
plot(x,dens(temp,x),type="l")
all.equal(matdens(temp,x),matdnorm(x, temp$mean,temp$sd))

temp2 = unimix(c(0.5,0.5),c(-2,2),c(1,1))
plot(x,dens(temp2,x),type="l")
