set.seed(22222)
source('ash.R')
source('jash.R')

# N: num of obs
# n: num of replicates
simdata = function(N, Nnull, altmean, altsd, betasd, n) {
  null = c(rep(1, Nnull), rep(0, N - Nnull))
  beta = c(rep(0, Nnull), rnorm(N - Nnull, altmean, altsd))
  Y = matrix(rnorm(N*n, rep(beta,n), rep(betasd,n)), ncol=n)
  return(list(null = null, beta = beta, Y = Y))
}

ss = simdata(10000, 8000, 0, 3, 1, 3)
beta.jash = jash(ss$Y)
beta.jash.auto = jash(ss$Y,auto=TRUE)

betahat = apply(ss$Y,1,mean)
sebetahat = apply(ss$Y,1,sd)/sqrt(dim(ss$Y)[2])
beta.ash.auto = ash(betahat,sebetahat,auto=TRUE)

# Histogram of true beta
hist(ss$beta, prob = TRUE, breaks = seq(-8, 8, length = 20))

plot(betahat,beta.jash$PosteriorMean,xlab='Observed betahat',ylab='Posterior mean',cex=0.5,pch=20)
points(betahat,beta.ash.auto$PosteriorMean, col=2, cex=0.5,pch=20)
title('Black: jash; Red: ash.auto')

# Compare sampled distribution P(beta|Y) to the normal mixture distribution P(beta|Y,hat(tau))
# Red lines are the conditional posterior distribution P(beta|Y,hat(tau)) 
par(mfrow=c(1,2))
sample.jash = posterior_sample_jash(beta.jash$post,1000,1:2)
hist(sample.jash$beta[1,],prob = TRUE, breaks = seq(-2, 2, length = 20),ylim=c(0,2.5),main='First unit')
x = seq(-2, 2, length = 1000)
normdens = dnorm(outer(x,rep(1,length(beta.jash$a.vec))),mean=beta.jash$condpost$normmean[rep(1,1000),],sd=beta.jash$condpost$normsd[rep(1,1000),])
mixnormdens = rowSums(beta.jash$condpost$condpi[rep(1,1000),]*normdens)
lines(x,mixnormdens,type='l',col=2) 

hist(sample.jash$beta[2,],prob = TRUE, breaks = seq(-2.5, 2.5, length = 20),ylim=c(0,2.5),main='Second unit')
x = seq(-2.5, 2.5, length = 1000)
normdens = dnorm(outer(x,rep(1,length(beta.jash$a.vec))),mean=beta.jash$condpost$normmean[rep(2,1000),],sd=beta.jash$condpost$normsd[rep(2,1000),])
mixnormdens = rowSums(beta.jash$condpost$condpi[rep(2,1000),]*normdens)
lines(x,mixnormdens,type='l',col=2) 

# Compare ash.auto, jash, jash.auto
o.ash = order(beta.ash.auto$qval)
o.jash = order(beta.jash$qval)
o.jash.auto = order(beta.jash.auto$qval)
par(mfrow=c(1,1))
plot(cumsum(ss$null[o.jash])/(1:10000),type='l',ylim=c(0,1),ylab='FDR' )
lines(cumsum(ss$null[o.ash])/(1:10000), col=2)
lines(cumsum(ss$null[o.jash.auto])/(1:10000), col='blue')
lines(sort(beta.jash$qvalue),lty=2)
lines(sort(beta.ash.auto$qvalue),lty=2,col=2)
lines(sort(beta.jash.auto$qvalue),lty=2,col='blue')
title('Dotted line: estimated FDR; Solid line: actual FDR')
legend('topleft',col=c('black','blue','red'),lty=c(1,1,1),legend=c('jash','jash.auto','ash.auto'))

beta.ash.auto$fitted
beta.jash$fitted
bata.jash.auto$fitted

# Use point mass
beta.jash.pm = jash(ss$Y,usePointMass=TRUE)
#beta.ash.pm = ash(betahat,sebetahat,auto=TRUE,usePointMass=TRUE)
beta.jash$fitted
beta.jash.pm$fitted
o.jash.pm = order(beta.ash.pm$qvalue)
plot(cumsum(ss$null[o.jash])/(1:10000),type='l',ylim=c(0,1),ylab='FDR' )
lines(cumsum(ss$null[o.jash.pm])/(1:10000), col=2)
lines(sort(beta.jash$qvalue),lty=2)
lines(sort(beta.jash.pm$qvalue),lty=2,col=2)
title('Dotted line: estimated FDR; Solid line: actual FDR')
legend('topright',col=c('black','red'),lty=c(1,1),legend=c('jash','jash.pm'))


# Use df modification for ash
beta.ash.t = ash(betahat,sebetahat,auto=TRUE,df=2)
o.ash.t = order(beta.ash.t$qval)
plot(cumsum(ss$null[o.jash])/(1:10000),type='l',ylim=c(0,1),ylab='FDR' )
lines(cumsum(ss$null[o.ash.t])/(1:10000), col=2)
lines(sort(beta.jash$qvalue),lty=2)
lines(sort(beta.ash.t$qvalue),lty=2,col=2)
title('Dotted line: estimated FDR; Solid line: actual FDR')
legend('topright',col=c('black','red'),lty=c(1,1),legend=c('jash','ash.t'))

# Large number of replicates
ss.largen = simdata(10000, 8000, 0, 3, 1, 20)
beta.largen.jash = jash(ss.largen$Y)
beta.largen.jash.auto = jash(ss.largen$Y,auto=TRUE)

betahat.largen = apply(ss.largen$Y,1,mean)
sebetahat.largen = apply(ss.largen$Y,1,sd)/sqrt(dim(ss.largen$Y)[2])
beta.largen.ash = ash(betahat.largen,sebetahat.largen,auto=TRUE)

o.largen.ash = order(beta.largen.ash$qval)
o.largen.jash = order(beta.largen.jash$qvalue)
o.largen.jash.auto = order(beta.largen.jash.auto$qvalue)
plot(cumsum(ss$null[o.largen.jash])/(1:10000),type='l',ylim=c(0,1),ylab='FDR' )
lines(cumsum(ss$null[o.largen.ash])/(1:10000), col=2)
lines(cumsum(ss$null[o.largen.jash.auto])/(1:10000), col='blue')
lines(sort(beta.largen.jash$qvalue),lty=2)
lines(sort(beta.largen.ash$qvalue),lty=2,col=2)
lines(sort(beta.largen.jash.auto$qvalue),lty=2,col='blue')
title('Large num of replicates case')
legend('topleft',col=c('black','blue','red'),lty=c(1,1,1),legend=c('jash','jash.auto','ash.auto'))


# Try some different data
truenull = c(rep(0, 1000), rep(1, 9000))
beta = c(rnorm(1000, -3, 1), rep(0, 9000))
s = rep(1, 10000)
Y = matrix(rnorm(30000, rep(beta,3), rep(s,3)),ncol=3)
Y.betahat = apply(Y,1,mean)
Y.sebetahat = apply(Y,1,sd)/sqrt(dim(Y)[2])
Y.ash.auto = ash(Y.betahat,Y.sebetahat,auto=TRUE)
Y.jash = jash(Y)
Y.jash.auto = jash(Y,auto=TRUE)

o = order(Y.jash.auto$qvalue)
plot(cumsum(truenull[o])/(1:10000), Y.jash.auto$qvalue[o], col = 2, type = "l",ylim=c(0,0.9))
#lines(cumsum(truenull[o])/(1:10000), Y.ash.auto$qval[o])
abline(a = 0, b = 1)

# Binomial data? UNDONE!
simdata.binom = function(N, Nnull, altmean, altsd, size, n) {
  null = c(rep(1, Nnull), rep(0, N - Nnull))
  beta = c(rep(0, Nnull), rnorm(N - Nnull, altmean, altsd))
  Y = matrix(rbinom(N*n, size=rep(size,N*n), prob=rep((beta+0.001)/size,n)), ncol=n)
  return(list(null = null, beta = beta, Y = Y))
}
bn = simdata.binom(10000, 8000, 0, 3, 5, 3)
