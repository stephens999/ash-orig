Model
========================================================
The observed values is in a $N\times n$ matrix $Y$ ($N$ tests and $n$ replicate observations for each tests), where \[y_{jk}|\beta_j,\tau_j\sim N(\beta_j,1/\tau_j),\ k=1,2,...,n,\ j=1,2,...,N\]
We assume that $\beta_j$ (mean), $\tau_j$ (precision) come from the Normal-Gamma mixture prior with $M\times L$ components:
\[\beta_j|\tau_j\sim N(\mu_j,1/(c_l\tau_j)),\quad \tau_j\sim Gamma(a_m,b_m)\quad w.p. \pi_{ml}\]
where $\sum_{m,l}\pi_{ml}=1,\ m=1,...,M,\ l=1,...,L$. Here $a_m>0,b_m>0,c_l>0,\mu_j$ are known constants (default setting $\mu_j=0$). Note that $c_l$ can be chosen as $\infty$, and thus result in a point mass at zero in $\beta$'s prior.

Posterior distribution
=======================
Hence, the conjugate postrior dis
tribution of $(\beta_j,\tau_j)$ given $Y$ ias also a Normal-Gamma mixture distribution, i.e. $\tau_j|Y\sim Gamma(\tilde{a},\tilde{b})$, $\beta_j|\tau_j,Y\sim N(\tilde{mu},1/\tilde{\tau})$. The full posterior distribution for $\beta_j$ is $P(\beta_j|Y)=\int P(\beta_j|\tau_j,Y)P(\tau_j|Y)d\tau_j$, which is hard to compute. For convenience, we can try using the distribution of $\beta_j|\hat{\tau_j},Y$ for further inference, where $\hat{\tau_j}=E(\tau_j|Y)$. Later we will compare the sampled $P(\beta_j|Y)$ to $P(\beta_j|\hat{\tau_j},Y)$.

Note: $P(\beta_j|Y)$ is feasible to compute, though very messy. I am working on it...

Estimate hyper-parameters
=======================
We estimate the hyper-parameters $\pi$ by maximizing the likelihood $L(\pi;Y)=\int P(Y,\beta,\tau|Y)d\beta d\tau$. 

Note: I find that sometimes if $a_m,b_m,c_l$ are not well chosen (e.g. the resulting $Var(Y|\beta,\tau)$ for each components overlap), then it is hard for the EM algorithm to converge. 

FDR and q-values
================
Same as ash, we use false sign rate for inference. 

Some examples
=============


```r
setwd("/Volumes/PERSONAL/MS")
set.seed(22222)
source("ash.R")
source("jash.R")

# simulate N observations with n replicates for each, Nnull under the null
# with altmean and altsd being the mean and sd of beta under the
# alternative
simdata = function(N, Nnull, altmean, altsd, betasd, n) {
    null = c(rep(1, Nnull), rep(0, N - Nnull))
    beta = c(rep(0, Nnull), rnorm(N - Nnull, altmean, altsd))
    Y = matrix(rnorm(N * n, rep(beta, n), rep(betasd, n)), ncol = n)
    return(list(null = null, beta = beta, Y = Y))
}
```


Simulate 10000 tests, with 2000 alternatives, β∼N(0,sd=3), and 8000 null values (β=0).

```r
ss = simdata(10000, 8000, 0, 3, 1, 3)
beta.jash = jash(ss$Y)
beta.jash.auto = jash(ss$Y, auto = TRUE)
```


We use ash for comparison. Here we use the test-wise sample-mean and sqrt(sample-variance/n) as the inputs of ash. (So the sebetahat here is not the true standard error of betahat!)

```r
betahat = apply(ss$Y, 1, mean)
sebetahat = apply(ss$Y, 1, sd)/sqrt(dim(ss$Y)[2])
beta.ash.auto = ash(betahat, sebetahat, auto = TRUE)
```


The histogram of true $\beta$:

```r
hist(ss$beta, prob = TRUE, breaks = seq(-8, 8, length = 20))
```

```
## Error: some 'x' not counted; maybe 'breaks' do not span range of 'x'
```


Both jash and ash give the Posterior means $E(\beta|Y)$, compare the shrinked means using two methods:

```r
plot(betahat, beta.jash$PosteriorMean, xlab = "Observed betahat", ylab = "Posterior mean", 
    cex = 0.5, pch = 20)
points(betahat, beta.ash.auto$PosteriorMean, col = 2, cex = 0.5, pch = 20)
title("Black: jash; Red: ash.auto")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


Compare the sampled full posterior distribution $P(\beta|Y)$ to the conditional posterior distribution $P(\beta|\hat{tau},Y)$, where $\hat{tau}=E(\tau|Y)$. Take the first two tests for illustration. It seems that in this case, the conditional distribution fits the full distribution well. (However, I am worried that it will make difference when $\tau$ has large scale or wide range. NEED TESTING!)


```r
# Compare sampled distribution P(beta|Y) to the normal mixture
# distribution P(beta|Y,hat(tau)) Red lines are the conditional posterior
# distribution P(beta|Y,hat(tau))
par(mfrow = c(1, 2))
sample.jash = posterior_sample_jash(beta.jash$post, 1000, 1:2)
hist(sample.jash$beta[1, ], prob = TRUE, breaks = seq(-2, 2, length = 20), ylim = c(0, 
    2.5), main = "First unit")
x = seq(-2, 2, length = 1000)
normdens = dnorm(outer(x, rep(1, length(beta.jash$a.vec))), mean = beta.jash$condpost$normmean[rep(1, 
    1000), ], sd = beta.jash$condpost$normsd[rep(1, 1000), ])
mixnormdens = rowSums(beta.jash$condpost$condpi[rep(1, 1000), ] * normdens)
lines(x, mixnormdens, type = "l", col = 2)

hist(sample.jash$beta[2, ], prob = TRUE, breaks = seq(-2.5, 2.5, length = 20), 
    ylim = c(0, 2.5), main = "Second unit")
x = seq(-2.5, 2.5, length = 1000)
normdens = dnorm(outer(x, rep(1, length(beta.jash$a.vec))), mean = beta.jash$condpost$normmean[rep(2, 
    1000), ], sd = beta.jash$condpost$normsd[rep(2, 1000), ])
mixnormdens = rowSums(beta.jash$condpost$condpi[rep(2, 1000), ] * normdens)
lines(x, mixnormdens, type = "l", col = 2)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


Compare the actual FDR and estimated FDR for ash, jash, jash.auto (a way of automatically selecting $a_m$,$b_m$,$c_l$). It seems that jash performs a bit better than ash (which makes sense, since the estimated sebetahat for ash is not accurate).

```r
# Compare ash.auto, jash, jash.auto
o.ash = order(beta.ash.auto$qval)
o.jash = order(beta.jash$qval)
o.jash.auto = order(beta.jash.auto$qval)
par(mfrow = c(1, 1))
plot(cumsum(ss$null[o.jash])/(1:10000), type = "l", ylim = c(0, 1), ylab = "FDR")
lines(cumsum(ss$null[o.ash])/(1:10000), col = 2)
lines(cumsum(ss$null[o.jash.auto])/(1:10000), col = "blue")
lines(sort(beta.jash$qvalue), lty = 2)
lines(sort(beta.ash.auto$qvalue), lty = 2, col = 2)
lines(sort(beta.jash.auto$qvalue), lty = 2, col = "blue")
title("Dotted line: estimated FDR; Solid line: actual FDR")
legend("topleft", col = c("black", "blue", "red"), lty = c(1, 1, 1), legend = c("jash", 
    "jash.auto", "ash.auto"))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


The fitted $\pi$ for jash: (sdvec gives the theoretical $SD(y_{jk}|a_m,b_m,c_l)$ for each component)

```r
beta.jash$fitted
```

```
## $pi
##  [1]  1.382e-43  1.391e-01  5.201e-01  6.636e-82  3.938e-45  1.502e-06
##  [7]  1.440e-01  2.093e-75  3.150e-77  6.903e-21  1.795e-01  1.064e-60
## [13] 4.866e-164  3.603e-97  1.720e-02  3.731e-43
## 
## $a.vec
##  [1] 5 4 3 2 5 4 3 2 5 4 3 2 5 4 3 2
## 
## $b.vec
##  [1] 50.0 20.0  1.0  0.1 50.0 20.0  1.0  0.1 50.0 20.0  1.0  0.1 50.0 20.0
## [15]  1.0  0.1
## 
## $c.vec
##  [1] 10.00 10.00 10.00 10.00  1.00  1.00  1.00  1.00  0.10  0.10  0.10
## [12]  0.10  0.01  0.01  0.01  0.01
## 
## $sdvec
##  [1]  0.06633  0.08944  0.20976  0.63561  0.11726  0.15811  0.37081
##  [8]  1.12361  0.60553  0.81650  1.91485  5.80230  2.34521  3.16228
## [15]  7.41620 22.47221
```

```r
bata.jash.auto$fitted
```

```
## Error: object 'bata.jash.auto' not found
```


We can also have a point mass at zero for $\beta$'s prior, by adding components with $c_l=\infty$. Compare the jash results with or without point mass:

```r
# Use point mass
beta.jash.pm = jash(ss$Y, usePointMass = TRUE)
# beta.ash.pm = ash(betahat,sebetahat,auto=TRUE,usePointMass=TRUE)
beta.jash$fitted
```

```
## $pi
##  [1]  1.382e-43  1.391e-01  5.201e-01  6.636e-82  3.938e-45  1.502e-06
##  [7]  1.440e-01  2.093e-75  3.150e-77  6.903e-21  1.795e-01  1.064e-60
## [13] 4.866e-164  3.603e-97  1.720e-02  3.731e-43
## 
## $a.vec
##  [1] 5 4 3 2 5 4 3 2 5 4 3 2 5 4 3 2
## 
## $b.vec
##  [1] 50.0 20.0  1.0  0.1 50.0 20.0  1.0  0.1 50.0 20.0  1.0  0.1 50.0 20.0
## [15]  1.0  0.1
## 
## $c.vec
##  [1] 10.00 10.00 10.00 10.00  1.00  1.00  1.00  1.00  0.10  0.10  0.10
## [12]  0.10  0.01  0.01  0.01  0.01
## 
## $sdvec
##  [1]  0.06633  0.08944  0.20976  0.63561  0.11726  0.15811  0.37081
##  [8]  1.12361  0.60553  0.81650  1.91485  5.80230  2.34521  3.16228
## [15]  7.41620 22.47221
```

```r
beta.jash.pm$fitted
```

```
## $pi
##  [1] 4.207e-222  1.497e-01  1.587e-04  0.000e+00 1.910e-225  1.255e-14
##  [7]  5.111e-01  0.000e+00 6.522e-207  1.299e-21  1.326e-01  0.000e+00
## [13]  0.000e+00  7.779e-77  1.903e-01 1.030e-297  0.000e+00  0.000e+00
## [19]  1.617e-02 1.841e-216
## 
## $a.vec
##  [1] 5 4 3 2 5 4 3 2 5 4 3 2 5 4 3 2 5 4 3 2
## 
## $b.vec
##  [1] 50.0 20.0  1.0  0.1 50.0 20.0  1.0  0.1 50.0 20.0  1.0  0.1 50.0 20.0
## [15]  1.0  0.1 50.0 20.0  1.0  0.1
## 
## $c.vec
##  [1]   Inf   Inf   Inf   Inf 10.00 10.00 10.00 10.00  1.00  1.00  1.00
## [12]  1.00  0.10  0.10  0.10  0.10  0.01  0.01  0.01  0.01
## 
## $sdvec
##  [1]  0.06325  0.06633  0.08944  0.20976  0.63561  0.11180  0.11726
##  [8]  0.15811  0.37081  1.12361  0.57735  0.60553  0.81650  1.91485
## [15]  5.80230  2.23607  2.34521  3.16228  7.41620 22.47221
```

```r
o.jash.pm = order(beta.ash.pm$qvalue)
```

```
## Error: object 'beta.ash.pm' not found
```

```r
plot(cumsum(ss$null[o.jash])/(1:10000), type = "l", ylim = c(0, 1), ylab = "FDR")
lines(cumsum(ss$null[o.jash.pm])/(1:10000), col = 2)
```

```
## Error: object 'o.jash.pm' not found
```

```r
lines(sort(beta.jash$qvalue), lty = 2)
lines(sort(beta.jash.pm$qvalue), lty = 2, col = 2)
title("Dotted line: estimated FDR; Solid line: actual FDR")
legend("topright", col = c("black", "red"), lty = c(1, 1), legend = c("jash", 
    "jash.pm"))
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


For ash, maybe the ``df'' option can resolve the understimating FDR issue? It works to some extent, but still not very satisfying.

```r
# Use df modification for ash
beta.ash.t = ash(betahat, sebetahat, auto = TRUE, df = 2)
o.ash.t = order(beta.ash.t$qval)
plot(cumsum(ss$null[o.jash])/(1:10000), type = "l", ylim = c(0, 1), ylab = "FDR")
lines(cumsum(ss$null[o.ash.t])/(1:10000), col = 2)
lines(sort(beta.jash$qvalue), lty = 2)
lines(sort(beta.ash.t$qvalue), lty = 2, col = 2)
title("Dotted line: estimated FDR; Solid line: actual FDR")
legend("topright", col = c("black", "red"), lty = c(1, 1), legend = c("jash", 
    "ash.t"))
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 


If the there are many replicates for each test, ash should yield pretty good results? Yes, it is. (But I am trying to figure out why jash performs this bad, maybe it has something to do with the settings of $a_m,b_m,c_l$)

```r
ss.largen = simdata(10000, 8000, 0, 3, 1, 20)
beta.largen.jash = jash(ss.largen$Y)
beta.largen.jash.auto = jash(ss.largen$Y, auto = TRUE)

betahat.largen = apply(ss.largen$Y, 1, mean)
sebetahat.largen = apply(ss.largen$Y, 1, sd)/sqrt(dim(ss.largen$Y)[2])
beta.largen.ash = ash(betahat.largen, sebetahat.largen, auto = TRUE)

o.largen.ash = order(beta.largen.ash$qval)
o.largen.jash = order(beta.largen.jash$qvalue)
o.largen.jash.auto = order(beta.largen.jash.auto$qvalue)
plot(cumsum(ss$null[o.largen.jash])/(1:10000), type = "l", ylim = c(0, 1), ylab = "FDR")
lines(cumsum(ss$null[o.largen.ash])/(1:10000), col = 2)
lines(cumsum(ss$null[o.largen.jash.auto])/(1:10000), col = "blue")
lines(sort(beta.largen.jash$qvalue), lty = 2)
lines(sort(beta.largen.ash$qvalue), lty = 2, col = 2)
lines(sort(beta.largen.jash.auto$qvalue), lty = 2, col = "blue")
title("Large num of replicates case")
legend("topleft", col = c("black", "blue", "red"), lty = c(1, 1, 1), legend = c("jash", 
    "jash.auto", "ash.auto"))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 


Try some different data: the issue of underestimating FDR is still noticable.

```r
truenull = c(rep(0, 1000), rep(1, 9000))
beta = c(rnorm(1000, -3, 1), rep(0, 9000))
s = rep(1, 10000)
Y = matrix(rnorm(30000, rep(beta, 3), rep(s, 3)), ncol = 3)
Y.betahat = apply(Y, 1, mean)
Y.sebetahat = apply(Y, 1, sd)/sqrt(dim(Y)[2])
Y.ash.auto = ash(Y.betahat, Y.sebetahat, auto = TRUE)
Y.jash = jash(Y)
Y.jash.auto = jash(Y, auto = TRUE)

o = order(Y.jash.auto$qvalue)
plot(cumsum(truenull[o])/(1:10000), Y.jash.auto$qvalue[o], col = 2, type = "l", 
    ylim = c(0, 0.9))
abline(a = 0, b = 1)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 


Works undone
=============
* Derive full posterior distribution $P(\beta|Y)$
* Better ways to (automatically) select $a_m,b_m,c_l$?
* Issue of underestimating false sign rate
* Use other approaches to estimate $\pi$ (e.g. VB?)
* Try more data (with different scales, underlying distributions, etc)
