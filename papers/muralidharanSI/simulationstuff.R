library(EbayesThresh)
library(locfdr)

"simulation1full" <- 
function(nreps = 100, twosided = FALSE)
{
	vecnorm = function(x,p=2) sum(x^p)^(1/p)
	# to avoid clutter		
			"simcomp1" <- 
			function(method, zmean, z, fullcomp = T, perrnorm = 2, ...)
			{
				nreps <- dim(z)[2]
				zs <- apply(z, 2, method, ...)
				mse <- apply(zs - zmean, 2, vecnorm, p = perrnorm)^perrnorm
				cat(method, "\n")
				return(mse)
			}
			
			"bayes" <- function(z, delta){
				N = length(z)
				e = list(pi = rep(1/N,N), mu = delta, sig = rep(1,N))
				return(effectSize(z,e)[,1])
			}
			
			"mixmodelfn" <-
			function(x,N,J){
				e = mixmodel(x,J=J, N = N, theonull = TRUE, sigInt = 1)
				return(effectSize(x,e)[,1])	
			}
			
			"mixmedian" = function(x,N,J){
				e = mixmodel(x,J=J, N = N, theonull = TRUE, sigInt = 1)
				return(postMedian(x,e))	
			}
			
			"locfdr2"<-
			function(z,df=5){
				bins = seq(min(z)-.1,max(z)+.1, len = 100)
				h = hist(z, bins, plot = F)
				x = h$m
				g  = glm(h$c~ns(x,df=df),family = poisson)
				ss = splinefun(x,log(dnorm(x)/g$fit),method = "natural")
				return(-ss(z,deriv=1))		
			}
			
			"fdrsmooth" <- 
			function(x, q = 0.050000000000000003)
			{
				#  implement the FDR procedure
				#
				n <- length(x)
				y <- rev(sort(abs(x)))
				t <-  - qnorm(((1:n) * q)/(2 * n))
				tfdr <- min(t[y >= t], t[1])
				xfdr <- x
				xfdr[abs(x) <= tfdr] <- 0
				return(xfdr)
			}

			"suresmooth" <- 
			function(x, adapt = T, returnthresh = F)
			{
				#
				#  find the SURE threshold and apply it to the sequence x
				#  also find carry out the ADAPT smoothing if adapt=T
				#
				n <- length(x)
				#
				#  first do the adapt test 
				#
				if(adapt & (sum(x^2)/n <= 1 + (log(n)/log(2))^1.5/sqrt(n))) 
						thresh <- sqrt(2 * log(n)) else {
					#
					#  the sure calculation
					#
					sx <- sort(x^2)
					sureseq <- n - 2 * (1:n) + seq(from = n - 1, to = 0,
						by = -1) * sx + cumsum(sx)
					thresh <- sqrt(min(sx[sureseq == min(sureseq)], 2 *
						log(n)))
				}
				xhat <- sign(x) * pmax(0, abs(x) - thresh)
				if(returnthresh)
					return(thresh)
				else return(xhat)
			}
			
			"universalsmooth" <- 
			function(x, soft = T)
			{
				n <- length(x)
				tt <- sqrt(2 * log(n))
				sx <- sign(x)
				ax <- abs(x)
				if(soft)
					xt <- sx * pmax(0, ax - tt)
				else {
					xt <- x
					xt[ax <= tt] <- 0
				}
				return(xt)
			}
	kvals <- c(5, 50, 500)
	muvals <- c(2, 3, 4, 5)
	results <- array(NA, dim = c(23, 4, 3, nreps), dimnames = list(
		c("exponential", "postmean", "SURE", 
		"FDR q=0.01", "FDR q=0.1", "FDR q=0.4", "universal soft", 
		"universal hard", "mixmodel J=3", "mixmodel N=10,J=3", "mixmodel N=50,J=3", "mixmodel N=100,J=3", "mixmodel N=200,J=3","mixmodel J=10", "mixmodel N=10,J=10","mixmodel N=50,J=10", "mixmodel N=100,J=10","mixmodel N=200,J=10", "locfdr7", "locfdr5", "mixmedian N=50,J=3", "mixmedian N=50,J=10","bayes"), as.character(muvals), as.character(
		kvals), NULL))
	#  simulate matrix of errors
	#
	set.seed(153)
	errmat <- matrix(rnorm(nreps * 1000), nrow = 1000)
	jperm <- sample(1:1000)
	#
	#  now conduct individual simulations
	#
	# if you want to calculate errors in a different norm, use the perrnorm parameter
	#  in simcomp1
	#
	for(jk in (1:3)) {
		for(jmu in (1:4)) {
			zmean <- rep(0, 1000)
			zmean[1:kvals[jk]] <- muvals[jmu]


			if(twosided == TRUE) zmean[1:floor(kvals[jk]/3)] = -zmean[1:floor(kvals[jk]/3)] 
					
			z <- matrix(rep(zmean, nreps), nrow = 1000) +
				errmat
			results[1, jmu, jk,  ] <- simcomp1(
				"ebayesthresh", zmean, z, a = NA, sd = 1)
			results[2, jmu, jk,  ] <- simcomp1(
				"ebayesthresh", zmean, z, threshrule = 
				"mean", a = NA, sd = 1)
			results[3, jmu, jk,  ] <- simcomp1("suresmooth",
				zmean, z, adapt = F)
			results[4, jmu, jk,  ] <- simcomp1("fdrsmooth",
				zmean, z, q = 0.01)
			results[5, jmu, jk,  ] <- simcomp1("fdrsmooth",
				zmean, z, q = 0.10000000000000001)
			results[6, jmu, jk,  ] <- simcomp1("fdrsmooth",
				zmean, z, q = 0.40000000000000002)
			results[7, jmu, jk,  ] <- simcomp1(
				"universalsmooth", zmean, z, soft = T)
			results[8, jmu, jk,  ] <- simcomp1(
				"universalsmooth", zmean, z, soft = F)
			results[9, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = -1, J = 3)
			results[10, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = 10, J = 3)
			results[11, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = 50, J = 3)
			results[12, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = 100, J = 3)
			results[13, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = 200, J = 3)
			results[14, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = -1, J = 10)
			results[15, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = 10, J = 10)
			results[16, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = 50, J = 10)
			results[17, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = 100, J = 10)
			results[18, jmu, jk,  ] <- simcomp1(
				"mixmodelfn", zmean, z, N = 200, J = 10)
			results[19, jmu, jk,  ] <- simcomp1(
				"locfdr2", zmean, z, df = 7)
			results[20, jmu, jk,  ] <- simcomp1(
				"locfdr2", zmean, z, df = 5)	
			results[21, jmu, jk,  ] <- simcomp1(
				"mixmedian", zmean, z, N = 50, J = 3)
			results[22, jmu, jk,  ] <- simcomp1(
				"mixmedian", zmean, z, N = 50, J = 10)
			results[23, jmu, jk,  ] <- simcomp1(
				"bayes", zmean, z, delta = zmean)
		}
	}
	return(results)
}

"table1ebayess" <- 
function(zs = ebayess.res, tex = F, avcal = T, relcal = F)
{
	#   takes results of simulation2() or simulation3() and
	#   produces Table 2.   Can either produce
	#   an array or else latex output
	#   If avcal, calculates means, if avcal=F calculates std errors
	#   If relcal=T does things relative to results for first method considered
	#
	if(relcal) {
		for(jj in (2:dim(zs)[1]))
			zs[jj,  ,  ,  ] <- zs[jj,  ,  ,  ] - zs[1,
				,  ,  ]
		zs[1,  ,  ,  ] <- 0
	}
	if(avcal) {
		zr <- apply(zs, c(1, 2, 3), mean)
	}
	else {
		zr <- sqrt(apply(zs, c(1, 2, 3), var)/dim(zs)[4])
	}
	if(!tex) {
		return(zr)
	}
	zt <- round(cbind(zr[,  , 2]))
	zt <- apply(zt, 1, paste, collapse = " & ")
	zt <- paste(names(zt), zt, sep = " & ", collapse = " \\  ")
	return(zt)
}

"table2ebayess" <- 
function(zs = ebayess.res, tex = F)
{
	#   takes results and
	#   produces inefficiency table.   Can either produce
	#   an array or else latex output
	#
	nmethods <- dim(zs)[1]
	methnames <- dimnames(zs)[[1]]
	zav <- apply(zs, c(1, 2, 3), mean)
	zav <- t(matrix(zav, nrow = nmethods, dimnames = list(methnames,
		NULL)))
	zbest <- apply(zav, 1, min)
	zineff <- (zav/zbest)
	zineffmat <- matrix(NA, nrow = nmethods, ncol = 5, dimnames = 
	list(methnames, c("mean", "median", "25th", "75th", "max")))
	zineffmat[, 1] <- apply(zineff, 2, mean)
	zineffmat[, 2] <- apply(zineff, 2, median)
	zineffmat[, 3] <- apply(zineff, 2, quantile, probs = 1/4)
	zineffmat[, 4] <- apply(zineff, 2, quantile, probs = 3/4)
	zineffmat[, 5] <- apply(zineff, 2, max)
	if(!tex) {
		return(zineffmat)
	}
	zt <- round(zineffmat)
	zt <- apply(zt, 1, paste, collapse = " & ")
	zt <- paste(names(zt), zt, sep = " & ", collapse = " \\  ")
	return(zt)
}