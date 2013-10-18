#setwd("~/Documents/git/ash/papers/muralidharanSI/")
library(locfdr)
library(fdrtool)
set.seed(355)

s = seq(0,5,by=0.01)
delta = c(rep(0,950), runif(50,2,4))
truemodel = list(pi = rep(1/1000,1000), mu = delta, sig = rep(1,1000))
nreps = 100
truefdr = FDRMixModel(s, truemodel, delta==0)$FDRT
errmat <- matrix(rnorm(nreps * 1000), nrow = 1000)
z = matrix(rep(delta,nreps),nrow=1000) + errmat
mix.fdr.fun = function(z,N,J,empnull=FALSE){
	theonull = !empnull
	a = NA
	if(theonull) a = 1
	e = mixFdr(z,J=J, P = N, theonull = theonull, noiseSD = a, calibrate = FALSE, plots = FALSE, nocheck = TRUE)
	ind = abs(e$mu - e$mu[1])<=.1
#	return(fdrMixModel(s,e,ind))
	return(FDRMixModel(s, e, ind)$FDRTwoSided)
}
logit = function(x) log(x/(1-x))
logistic = function(x) 1/(1+exp(-x))
loc.fdr.fun = function(z,df=7,nulltype=0){
#	l = locfdr(z, df = df, nulltype = nulltype, plot = 0)
#	find = 8 * (nulltype == 0) + 2 * (nulltype != 0)
#	sp = splinefun(l$mat[,1], l$mat[,find], method = "natural")
#	pmin(1,pmax(0, sp(s)))
	l = locfdr(z, df = df, nulltype = nulltype, plot = 0, sw = 2)
	sp = splinefun(l$x, logit(pmin(cumsum(l$f+0.001)/sum(l$f+0.001),0.999)), method = "natural")
	pmin(1,pmax(0, l$pds[1]*((s<0) + (s>=0)*(1 - pnorm(s, l$pds[2], l$pds[3]) + pnorm(-s,l$pds[2],l$pds[3])) / (1 - logistic(sp(s)) + logistic(sp(-s))))))
}
fdrtool.fdr.fun = function(z, empnull = F){
	f = fdrtool(z, "normal",plot=FALSE)	
	res = f$lfdr		
	if(empnull == F){
	f = fdrtool(2*(1 - pnorm(abs(z))), "pvalue", plot = FALSE)
	res = f$qval	
	}
	sp = splinefun(z, logit(pmax(pmin(res,0.999),.001)), method = "natural")	
	pmin(1,pmax(0,logistic(sp(s))))
}

wt = dens.mixture(s, truemodel)

resMfdrTheo = apply(z, 2, mix.fdr.fun, N = 50, J = 10)
resLocfdrTheo = apply(z, 2, loc.fdr.fun)
resFdrtoolTheo = apply(z, 2, fdrtool.fdr.fun)

resMfdrEmp = apply(z, 2, mix.fdr.fun, N = 50, J = 10, empnull = TRUE)
resLocfdrEmpMLE = apply(z, 2, loc.fdr.fun, nulltype = 1)
resLocfdrEmpGeo = apply(z, 2, loc.fdr.fun, nulltype = 2)
resFdrtoolEmp = apply(z, 2, fdrtool.fdr.fun, empnull = TRUE)

biasCurveMfdrTheo = apply((resMfdrTheo), 1, mean)
biasCurveLocfdrTheo = apply((resLocfdrTheo+0.00001), 1, mean)
biasCurveFdrtoolTheo = apply((resFdrtoolTheo+0.00001), 1, mean)
biasCurveMfdrEmp = apply((resMfdrEmp), 1, mean)
biasCurveLocfdrEmpM = apply((resLocfdrEmpMLE), 1, mean)
biasCurveLocfdrEmpG = apply((resLocfdrEmpGeo), 1, mean)
biasCurveFdrtoolEmp = apply((resFdrtoolEmp), 1, mean)

op = par(mfrow = c(1,2), xpd = FALSE, mar = c(5,2,2,2)+0.1)
plot(s, (truefdr), t = 'l', col = "black", lwd = 2, ylim = c(0,1), xlab = "z", ylab = "E(FDR.hat(z))", main = "Expected FDR estimates")
lines(s, (biasCurveFdrtoolEmp), col = 4, lty = 2, lwd = 1.5)
lines(s, (biasCurveFdrtoolTheo), col = 4, lwd = 1.5)
lines(s, (biasCurveLocfdrEmpM), col = 3, lty = 2, lwd = 1.5)
lines(s, (biasCurveLocfdrEmpG), col = 5, lty = 2, lwd = 1.5)
lines(s, (biasCurveLocfdrTheo), col = 3, lwd = 1.5)
lines(s, (biasCurveMfdrTheo), col = 2, t='l', lwd = 1.5)
lines(s, (biasCurveMfdrEmp), col = 2, t='l', lty = 2, lwd = 1.5)
lines(s, truefdr, col = "black", lwd = 2)
names = c("MixFdr Th", "Locfdr Th", "Fdrtool Th", "MixFdr Emp", "Locfdr MLE", "Locfdr CM", "Fdrtool Emp")
lcols = c(2,3,4,2,3,5,4)
ltps = c(1,1,1,2,2,2,2)
#names = c("MixModel Th", "Locfdr Th", "Mixmodel Emp", "Locfdr MLE")
#lcols = c(2,3,2,3)
#ltps = c(1,1,2,2)

varLCurveMfdrTheo = apply((resMfdrTheo), 1, sd)
varLCurveLocfdrTheo = apply((resLocfdrTheo), 1, sd)
varLCurveFdrtoolTheo = apply((resFdrtoolTheo), 1, sd)
varLCurveMfdrEmp = apply((resMfdrEmp), 1, sd)
varLCurveLocfdrEmpM = apply((resLocfdrEmpMLE), 1, sd)
varLCurveLocfdrEmpG = apply((resLocfdrEmpGeo), 1, sd)
varLCurveFdrtoolEmp = apply((resFdrtoolEmp), 1, sd)
par(mar = c(5,2,2,8)+0.1)
plot(s, varLCurveMfdrTheo, col = 2, t='l', xlab = "z", ylab = "Sd(FDR.hat(z))",main = "Standard Deviation of FDR estimates", lwd = 1.5, ylim = c(0,0.1))
lines(s, varLCurveFdrtoolTheo, col = 4, ylim=c(0,2), lwd = 1.5)
lines(s, varLCurveMfdrEmp, col = 2, t='l',lty=2, lwd = 1.5)
lines(s, varLCurveFdrtoolEmp, col = 4, ylim=c(0,2),lty=2, lwd = 1.5)
lines(s, varLCurveLocfdrTheo, col = 3, ylim=c(0,2), lwd = 1.5)
lines(s, varLCurveLocfdrEmpM, col = 3, ylim=c(0,2),lty=2, lwd = 1.5)
lines(s, varLCurveLocfdrEmpG, col = 5, ylim=c(0,2),lty=2, lwd = 1.5)
par(xpd = TRUE)
legend(5.1, .105, legend = c("Truth",names), col = c("black",lcols), lty = c(1,ltps), lwd = c(2, rep(1.5,length(names))), bty='n')
par(op)


# see how well they estimate the threshholds
thresh = function(fdr){
	q = seq(0.01,0.20,by=0.01)
	o = outer(fdr, q, FUN = ">")
	inds = apply(apply(o, 2, cumsum),2,which.max)
	s[inds]
}

q = seq(0.01,0.20,by=0.01)
tMfdrTh = apply(resMfdrTheo, 2, thresh)
tLocfdrTh = apply(resLocfdrTheo, 2, thresh)
tFdrtoolTh = apply(resFdrtoolTheo, 2, thresh)

tMfdrEmp = apply(resMfdrEmp, 2, thresh)
tLocfdrEmpMLE = apply(resLocfdrEmpMLE, 2, thresh)
tLocfdrEmpGeo = apply(resLocfdrEmpGeo, 2, thresh)
tFdrtoolEmp = apply(resFdrtoolEmp, 2, thresh)

op = par(mfrow = c(1,2), xpd = FALSE, mar = c(5,2,2,2)+0.1)
truethresh = thresh(truefdr)
plot(q, truethresh, t = 'l', lwd =2, ylab = "E(t.hat(q))", main = "Expected Threshhold Estimate", ylim = c(2.5,4.5))
lines(q,apply(tMfdrTh,1,mean),col=2, lwd = 1.5)
lines(q,apply(tLocfdrTh,1,mean),col=3, lwd = 1.5)
lines(q,apply(tFdrtoolTh,1,mean),col=4, lwd = 1.5)
lines(q,apply(tMfdrEmp,1,mean),col=2, lty = 2, lwd = 1.5)
lines(q,apply(tLocfdrEmpMLE,1,mean),col=3, lty = 2, lwd = 1.5)
lines(q,apply(tLocfdrEmpGeo,1,mean),col=5, lty = 2, lwd = 1.5)
lines(q,apply(tFdrtoolEmp,1,mean),col=4, lty = 2, lwd = 1.5)
par(xpd = TRUE)
par(xpd = FALSE, mar = c(5,2,2,8)+.1)
plot(q, apply(tMfdrTh,1,sd), col = 2, t = 'l', ylim = c(0,.5), lwd = 1.5, ylab = "Sd(t.hat(q))", main = "Standard Deviation of Threshold Estimate")
lines(q, apply(tLocfdrTh,1,sd), col = 3, t = 'l', lwd = 1.5)
lines(q, apply(tFdrtoolTh,1,sd), col = 4, t = 'l', lwd = 1.5)
lines(q, apply(tMfdrEmp,1,sd), col = 2, t = 'l', lty = 2, lwd = 1.5)
lines(q, apply(tLocfdrEmpMLE,1,sd), col = 3, t = 'l', lty = 2, lwd = 1.5)
lines(q, apply(tLocfdrEmpGeo,1,sd), col = 5, t = 'l', lty = 2, lwd = 1.5)
lines(q, apply(tFdrtoolEmp,1,sd), col = 4, t = 'l', lty = 2, lwd = 1.5)
par(xpd = TRUE)
legend(0.21,0.525, legend = c("Truth",names), col = c("black",lcols), lty = c(1,ltps), lwd = c(2, rep(1.5,length(names))), bty = 'n')
par(op)