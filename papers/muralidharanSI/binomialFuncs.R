binomgroupProbs = function(z, NN, e){
	phiMat = binomdensMat(z,NN,e)
	pi = e$pi
	dens = (phiMat %*% pi)[,1]
	gpProbs = phiMat / matrix(dens, nrow(phiMat), ncol(phiMat), byrow = F) * matrix(pi, nrow(phiMat), ncol(phiMat), byrow =T)
}

binomdensMat = function(z, NN, e){
	pi = e$pi
	betaparams = rbind(e$alpha,e$beta)
	J = length(pi)
	N = length(z)
	alphas = betaparams[1,]
	betas = betaparams[2,]
	alphas = matrix(alphas, N, J, T)
	betas = matrix(betas, N, J, T)
	Z = matrix(z,N,J)
	NNz = matrix(NN,N,J)
	densMat = lgamma(alphas + betas) - lgamma(alphas) - lgamma(betas)
	densMat = densMat + lgamma(Z + alphas) + lgamma(NNz - Z + betas) - lgamma(NNz + alphas + betas)
	densMat = densMat + lchoose(NNz,Z)
	densMat = exp(densMat)
}



emfitter = function(z, NN, J = 2, starts = NA, tol = 0.001, maxiter = 100){
	
	#put some functions in here so that they don't clutter up
	#the workspace
	# MLE is too hard to get so let's use MOM estimates

betaDataMOM = function(z, NN, Ze){
	rsums = colSums(Ze)	
	deltasur = z/NN
	xbars = (t(Ze) %*% deltasur)[,1]/rsums
	xsqbars = (t(Ze) %*% deltasur^2)[,1]/rsums
	ssqs = xsqbars - xbars^2
	muHats = xbars
	a1 = (t(Ze) %*% (1-1/NN))[,1]/rsums
	a2 = (t(Ze) %*% (1/NN))[,1]/rsums
	mu1min = muHats*(1-muHats)
	Mhats = 1 + mu1min * a1
	denom = pmax(0.001,ssqs - mu1min * a2)
	Mhats = Mhats/denom
	alphas = muHats * Mhats
	betas = (1-muHats) * Mhats
	return(rbind(alphas,betas))	
}

binomestep = function(z, pi, betaparams, NN){
	G = binomgroupProbs(z, NN, e = list(pi = pi, alpha = betaparams[1,], beta = betaparams[2,]))	
	N = length(z)
	J = length(pi)
	Z = matrix(z, length(z), length(pi))
	NNz = matrix(NN, length(z), length(pi))
	alphas = matrix(betaparams[1,], nrow = N, ncol = J, byrow = T)
	betas = matrix(betaparams[2,], nrow = N, ncol = J, byrow = T)
	deltaMat = (alphas + Z)/(alphas+betas+NNz)
	deltas = rowSums(deltaMat * G)		
	return(list(Ze = G, deltas = deltas))
}

	
	
	if(is.na(starts)){
		pi = rep(1/J,J)	
		alphas = rnorm(J, mean = 3, sd = 1)
		betas = 3 - alphas + rnorm(J, mean = 1, sd = 0.1)
		starts = cbind(pi, abs(alphas)+1,abs(betas)+1)
	}
	pi = starts[,1]
	betaparams = t(starts[,2:3])
	# that's it, do the loop	
	iter = 0
	converged = FALSE
	while(!converged && iter<maxiter){
		Zd = binomestep(z, pi, betaparams, NN)
		respi = colMeans(Zd$Ze)
		resbetaparams = betaDataMOM(z,NN,Zd$Ze) 
		distance = sum((pi-respi)^2)+sum((resbetaparams-betaparams)^2) 
		if(distance<tol) converged = TRUE
		iter = iter + 1		
		pi = respi
		betaparams = resbetaparams
	}
	return(list(pi = pi, alphas = betaparams[1,], betas = betaparams[2,], converged = converged, nIter = iter, z = z, NN = NN))
}


binomialEffectSize = function(z,e,NN){
	G = binomgroupProbs(z, NN, e)	
	N = length(z)
	J = length(e$pi)
	Z = matrix(z, N, J)
	NNz = matrix(NN, N, J)
	alphas = matrix(e$alpha, nrow = N, ncol = J, byrow = T)
	betas = matrix(e$beta, nrow = N, ncol = J, byrow = T)
	deltaMat = (alphas + Z)/(alphas+betas+NNz)
	deltas = rowSums(deltaMat * G)		
}

binomMixDens = function(i, e, NN){
	d = binomdensMat(i, NN, e)
	(d %*% e$pi)[,1]
}

EmuZ = function(z,NN,e){
	
	EmuBeta = function(ab){
		if(ab[1]/ab[2] < 0.001){
			return(0)
		}
		else{
			integrate(f,0,1,a=ab[1],b=ab[2])$value
		}	
	}
	f = function(x, a, b){
		asin(sqrt(x)) * dbeta(x,a,b)	
	}
	
	pi = e$pi
	alpha = e$alpha
	beta = e$beta
	J = length(pi)
	N = length(z)
	postAlphas = matrix(alpha, nrow = N, ncol = J, byrow = T)
	postAlphas = postAlphas + matrix(z, N, J)
	postBetas = matrix(beta, N,J,T)
	postBetas = postBetas + matrix(NN-z,N,J)
	postBetaParams = array(0,c(N,J,2))
	postBetaParams[,,1] = postAlphas
	postBetaParams[,,2] = postBetas
	postMeanMat = apply(postBetaParams,c(1,2),EmuBeta)
	G = binomgroupProbs(z, NN, e)
	postMeans = postMeanMat * G
	return(rowSums(postMeans))
}
