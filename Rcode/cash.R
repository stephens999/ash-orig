library(ordinal)
source('Rcode/ash.R')

# EXPECTATION STEP
expectation=function(X,likely,theta,beta){
	n=dim(X)[1]
	linear=matrix(rep(theta,n),n,byrow=T)-colSums(t(X)*beta)
	priorp=t(apply(cbind(rep(0,n),exp(linear)/(1+exp(linear)),rep(1,n)),1,diff))
	margin=sum(log(rowSums(likely*priorp)))
	poster=likely*priorp/rowSums(likely*priorp)
	return(list(poster=poster,margin=margin))}

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
	print(beta)
	return(list(theta=theta,beta=beta))}
	
# CASH
# GIVEN: ESTIMATED EFFECTS (alpha), STANDARD ERRORS (S), COVARIATE MATRIX (X) AND METHOD
# alpha: n-vector
# s: n-vector
# X: n-by-p matrix
# glm.method: 'clm' (default), 'plum' or 'vglm'

cash=function(alpha,s,X,ltol=.0001,maxiter=5000){
	
	n=dim(X)[1]
	p=dim(X)[2]
	
	# SPECIFY: VALUES OF SIGMA TO CONSIDER
	sigma=c(.001,.002,.004,.008,.016,.032,.064,.128,.256,.512,1.024,2.048,4.096) # K-vector
	
	# INITIALIZE PARAMETERS (INITIALIZE BETA TO BE ZERO AND FIT THETA USING ASH)
	beta=rep(0,p)# p-vector
	probs=ash(alpha,s,prior='uniform',sigmaavec=sigma)$fit$pi
	gamma=diffinv((probs+.01)/sum(probs+.01))[2:length(probs)]
	theta=log(gamma/(1-gamma))
#	theta=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
	maximized=list(theta=theta,beta=beta)
	
	varmat=2*(matrix(rep(sigma^2,n),n,byrow=T)+s^2)
	likely=exp(-alpha^2/varmat)/sqrt(pi*varmat)

	beta=rep(0,maxiter)
	theta=matrix(0,maxiter,length(sigma)-1)
	margin=rep(0,maxiter)
	l=0

	for(i in 1:maxiter){
		expected=expectation(X,likely,maximized$theta,maximized$beta)
		maximized=maximization(expected$poster,X,maximized$theta,maximized$beta)
		beta[i]=maximized$beta
		theta[i,]=maximized$theta
		margin[i]=expected$margin
		if(abs(expected$margin-l)<ltol)break;
		l=expected$margin
	}
	
	beta=beta[1:i]
	theta=theta[1:i,]
	margin=margin[1:i]
	
	return(list(maximized=maximized,beta=beta,theta=theta,margin=margin))
}