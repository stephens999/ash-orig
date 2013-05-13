library(ordinal)
library(VGAM)
library(catdata)
source('plum.R')
source('../ash.R')

# EXPECTATION STEP
expectation=function(X,likely,theta,alpha){
	n=length(beta)
	linear=matrix(rep(theta,n),n,byrow=T)-colSums(t(X)*alpha)
	priorp=t(apply(cbind(rep(0,n),exp(linear)/(1+exp(linear)),rep(1,n)),1,diff))
	margin=sum(log(colSums(likely*priorp)))
	print(margin)
	poster=likely*priorp/rowSums(likely*priorp)
	return(poster)}

# MAXIMIZATION STEP (CLM)
maximization.clm=function(poster,X){
	K=dim(poster)[2]
	n=dim(X)[1]
	p=dim(X)[2]
	output=ordered(rep(1:K,each=n))
	weight=array(poster)
	design=kronecker(array(1,K),X)
	fitted=clm(output~design,weights=weight)$coe
	theta=fitted[1:(K-1)]
	alpha=fitted[K:(K+p-1)]
	return(list(theta=theta,alpha=alpha))}

# MAXIMIZATION STEP (PLUM)
maximization.plum=function(poster,X){
	K=dim(poster)[2]
	p=dim(X)[2]
	fitted=plum(poster,X)$coef
	theta=fitted[1:(K-1)]
	alpha=-fitted[K:(K+p-1)]
	return(list(theta=theta,alpha=alpha))}
	
# MAXIMIZATION STEP (POLR) ?
	
# MAXIMIZATION STEP (VGLM)
maximization.vglm=function(poster,X){
	K=dim(poster)[2]
	n=dim(X)[1]
	p=dim(X)[2]
	fitted=Coef.vlm(vglm(round(1e16*poster)~X,family=propodds(F)))
	theta=fitted[1:(K-1)]
	alpha=-fitted[K:(K+p-1)]
	return(list(theta=theta,alpha=alpha))}
	
# CASH
# GIVEN: ESTIMATED EFFECTS (BETA), STANDARD ERRORS (S), COVARIATE MATRIX (X) AND METHOD
# beta: n-vector
# s: n-vector
# X: n-by-p matrix
# glm.method: 'clm' (default), 'plum' or 'vglm'

cash=function(beta,s,X,glm.method='clm'){
	
	n=length(beta)
	
	if(glm.method=='plum'){maximization=maximization.plum}
	# (n,p,K)=(100,2,13): 275 seconds
	if(glm.method=='clm'){maximization=maximization.clm}
	# (n,p,K)=(100,2,13): .64 seconds
	if(glm.method=='vglm'){maximization=maximization.vglm}
	# haven't gotten this one working yet

	# SPECIFY: VALUES OF SIGMA TO CONSIDER
#	sigma=c(.001,.002,.004,.008,.016,.032,.064,.128,.256,.512,1.024,2.048,4.096) # K-vector
	sigma= c(0.00025,0.0005,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1.024,2.048,4.096,8.192)
	
	# INITIALIZE PARAMETERS (INITIALIZE ALPHA TO BE ZERO AND FIT THETA USING ASH)
	alpha=0 # p-vector
	probs=ash(beta,s)$fit$pi
	gamma=diffinv((probs+.01)/sum(probs+.01))[2:length(probs)]
	theta=log(gamma/(1-gamma))
	theta=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
	maximized=list(theta=theta,alpha=alpha)
	
	varmat=2*(matrix(rep(sigma^2,n),n,byrow=T)+s^2)
	likely=exp(-beta^2/varmat)/sqrt(pi*varmat)

	for(i in 1:100){
		expected=expectation(X,likely,maximized$theta,maximized$alpha)
		maximized=maximization(expected,X)
	}
	return(maximized)
}

# SIMULATION
set.seed(4)
n=1000
X=cbind(matrix(c(rep(1,.1*n),rep(0,.9*n))),matrix(c(rep(0,.9*n),rep(1,.1*n))))
sigma=c(.001,.002,.004,.008,.016,.032,.064,.128,.256,.512,1.024,2.048,4.096)
theta=c(-3,-2.5,-2,-1.5,-1,-.5,.5,1,1.5,2,2.5,3)
alpha=c(-3,3)
s=exp(rnorm(n))
linear=matrix(rep(theta,n),n,byrow=T)-colSums(t(X)*alpha)
priorp=t(apply(cbind(rep(0,n),exp(linear)/(1+exp(linear)),rep(1,n)),1,diff))
simz=function(probs){rmultinom(1,1,probs)}
z=t(apply(priorp,1,simz))
sigmai=rowSums(z*matrix(rep(sigma,n),n,byrow=T))
beta=rnorm(n,0,sqrt(s^2+sigmai^2))

start=Sys.time()
cash(beta,s,X,glm.method='clm')
Sys.time()-start