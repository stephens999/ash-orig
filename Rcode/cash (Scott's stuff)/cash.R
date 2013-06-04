library(ordinal)
source('../ash.R')

# EXPECTATION STEP
expectation=function(X,likely,theta,beta){
	n=dim(X)[1]
	linear=matrix(rep(theta,n),n,byrow=T)-colSums(t(X)*beta)
	priorp=t(apply(cbind(rep(0,n),exp(linear)/(1+exp(linear)),rep(1,n)),1,diff))
	margin=sum(log(rowSums(likely*priorp)))
	print(c(margin, sum(priorp)))
	poster=likely*priorp/rowSums(likely*priorp)
	return(poster)}

# MAXIMIZATION STEP
maximization=function(poster,X,thetastart,betastart){
	K=dim(poster)[2]
	n=dim(X)[1]
	p=dim(X)[2]
	output=ordered(rep(1:K,each=n))
	weight=array(poster)
	design=kronecker(array(1,K),X)
	fitted=clm(output~design,weights=weight,start=c(thetastart,betastart))$coe
	theta=fitted[1:(K-1)]
	beta=fitted[K:(K+p-1)]
	return(list(theta=theta,beta=beta))}
	
# CASH
# GIVEN: ESTIMATED EFFECTS (alpha), STANDARD ERRORS (S), COVARIATE MATRIX (X) AND METHOD
# alpha: n-vector
# s: n-vector
# X: n-by-p matrix
# glm.method: 'clm' (default), 'plum' or 'vglm'

cash=function(alpha,s,X){
	
	n=dim(X)[1]
	p=dim(X)[2]
	
	# SPECIFY: VALUES OF SIGMA TO CONSIDER
#	sigma=c(.001,.002,.004,.008,.016,.032,.064,.128,.256,.512,1.024,2.048,4.096) # K-vector
	sigma= c(0.00025,0.0005,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1.024,2.048,4.096,8.192)
	
	# INITIALIZE PARAMETERS (INITIALIZE BETA TO BE ZERO AND FIT THETA USING ASH)
	beta=rep(0,p)# p-vector
	probs=ash(alpha,s)$fit$pi
	gamma=diffinv((probs+.01)/sum(probs+.01))[2:length(probs)]
	theta=log(gamma/(1-gamma))
#	theta=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
	maximized=list(theta=theta,beta=beta)
	
	varmat=2*(matrix(rep(sigma^2,n),n,byrow=T)+s^2)
	likely=exp(-alpha^2/varmat)/sqrt(pi*varmat)

	for(i in 1:100){
		expected=expectation(X,likely,maximized$theta,maximized$beta)
		maximized=maximization(expected,X,maximized$theta,maximized$beta)
	}
	return(maximized)
}