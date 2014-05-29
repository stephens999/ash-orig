# When measurement precision is constant, differences in $q$-values between methods should be driven by estimates of $\pi_0$

Here we explore the root cause of differences between $q$-values
produced by different methods, in the simple case where measurement precisions $(se(\hat\beta_j))$
are the same for all observations. We claim that these differences
should be driven primarily by the fact that different methods produce 
different estimates for $\pi_0$, rather
than other differences in methodology. 

First, we outline the intuition behind why differences should be driven primarily by differences in $\pi_0$. Suppose we are given a set of $n$ $p$-values corresponding to $n$ distinct null hypotheses.
Given an estimate of the proportion of true nulls, $\pi_0$, the false discovery
rate at threshold $p$, FDR($p$), is easily estimated "non-parametrically" using the fact that the expected number of false discoveries at threshold $p$ is 
$n p \pi_0$. The $q$-value for observation $j$ is then, by definition, $q_j = \text{FDR}(p_j)$.
This simple "non-parametric" approach to computing $q$-values is implemented in the following function:

```r
qval.nonp = function(p, pi0) {
    n = length(p)
    efp = n * p * pi0  # expected number of false positives at threshold p
    np = rank(p)  #number of positives
    return(efp/np)
}
```


This approach is sufficiently straightforward and natural that one might make the following claim:
any method that, using the same value $\pi_0$, produces appreciably different $q$-values
from the non-parametric approach, should be considered with suspicion (or at least with caution).
This would in turn imply that differences between (sensible) methods should be driven by differences in $\pi_0$, since if they have the same $\pi_0$ then they should both produce results similar to the simple nonparametric approach.

This said, the method for computing $q$-values in `ash` is, on the face of it, very different from the non-parametric approach, and so it is not clear in advance that the two will match closely. In `ash` $q$-values are computed by first computing the local false discovery rate (lfdr) for each observation, and then, at any given threhold, 
averaging the lfdr values that exceed this threshold to estimate the $q$-value at that threshold. 
This approach is effectively estimating a tail-area (the $q$ value) by first estimating
a density (the lfdr) and then approximating an integral (by averaging). 
From this viewpoint, `ash` tackles an easy problem (estimating a tail area) by first solving a harder problem (estimating a density), and one might worry that this could introduce
errors.

We examine this issue by simulation. First we define a simple simulation function that simulates
a mixture of 0 and non-zero (normally-distributed) beta values.

```r
# simulate n beta-hat values, nnull under the null with altmean and altsd
# being the mean and sd of beta under the alternative
simdata = function(n, nnull, altmean, altsd, betahatsd) {
    null = c(rep(1, nnull), rep(0, n - nnull))
    beta = c(rep(0, nnull), rnorm(n - nnull, altmean, altsd))
    betahat = rnorm(n, beta, betahatsd)
    return(list(null = null, beta = beta, betahat = betahat, betahatsd = betahatsd))
}
library(ashr)
```

```
## Loading required package: truncnorm
## Loading required package: SQUAREM
```

```r
library(qvalue)
library(locfdr)
```

```
## Loading required package: splines
```

```r
library(mixfdr)
```


In our first simulation, half the $\beta_j$ are null, and the other half are $\beta_j \sim N(0,2^2)$. Then $\hat\beta_j \sim N(\beta_j, 1)$. We apply `qvalue`,`ash`,`locfdr` and `mixFdr` to estimate the $q$-values. (The last three methods all estimate local fdr values; we convert these to corresponding estimated $q$ values using the `ashr` function `qval.from.lfdr`, which, for each observation $j$, computes the $q$ value as the average of the lfdr values more significant than observation $j$.)

```r
set.seed(100)
test = simdata(10000, 5000, 0, 2, 1)
test.ash = ash(test$betahat, test$betahatsd, method = "fdr")
z = test$betahat/test$betahatsd
p = pchisq(z^2, df = 1, lower.tail = FALSE)
test.qv = qvalue(p)
test.locfdr = locfdr(z, nulltype = 0, plot = 0)
test.locfdr.qval = qval.from.lfdr(test.locfdr$fdr)
test.mixfdr = mixFdr(z, theonull = TRUE, noiseSD = 1, plots = FALSE)
```

```
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Null proportion pi0 is small. Consider increasing penalization and/or using an empirical null.
## Warning: Using an empirical null with a fitted noiseSD gives a substantially different model. Consider rerunning with theonull = FALSE and noiseSD = NA.
```

```r
test.mixfdr.qval = qval.from.lfdr(test.mixfdr$fdr)
```


The estimates of $\pi_0$ from the four methods are given in the following table:

```r
kable(data.frame(method = c("qvalue", "ash", "locfdr", "mixfdr"), pi0 = round(c(test.qv$pi0, 
    get_pi0(test.ash), test.locfdr$fp0[1, 3], test.mixfdr$pi[1]), digits = 2)))
```

```
## |method  |   pi0|
## |:-------|-----:|
## |qvalue  |  0.72|
## |ash     |  0.47|
## |locfdr  |  0.73|
## |mixfdr  |  0.85|
```


Now we perform a "sanity check" for each method, comparing its $q$-values against the
nonparametric $q$ values obtained using the same estimate of $\pi_0$. Both `ash` and `qvalue`
clearly pass this sanity check (for `qvalue`, this is almost tautological, since it effectively uses the non-parametric procedure). For `locfdr`, it also passes the sanity check except for
some artifacts near $p=1$ that we do not fully understand. For `mixfdr` there seems to be a systematic deviation from the non-parametric method that suggests...

```r
par(mfcol = c(2, 2))
plot(qval.nonp(p, get_pi0(test.ash)), test.ash$qvalue, xlab = "nonparametric q-values (with pi0 estimated by ash)", 
    ylab = "ash q-values", main = "ash")
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(qval.nonp(p, test.qv$pi0), test.qv$qvalues, xlab = "nonparametric q-values (with pi0 estimated by qvalue)", 
    ylab = "qvalue q-values", main = "qvalue")
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(qval.nonp(p, test.mixfdr$pi[1]), test.mixfdr.qval, xlab = "nonparametric q-values (with pi0 estimated by mixfdr)", 
    ylab = "mixfdr q-values", main = "mixfdr")
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(qval.nonp(p, test.locfdr$fp0[1, 3]), test.locfdr.qval, xlab = "nonparametric q-values (with pi0 estimated by locfdr)", 
    ylab = "locfdr q-values", main = "locfdr")
abline(a = 0, b = 1, col = 2, lwd = 2)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

```r
par(mfcol = c(1, 1))
```




In this second simulation, half the $\beta_j$ are null, and the alternatives
are $\beta_j \sim N(2,1)$. Note that this deviates quite a bit from the `ash` assumption that
the non-zero $\beta_j$ will have a unimodal distribution centered on 0. This presumably
accounts for the slightly worse correspondance between the `ash` and non-parametric $q$-value compared with the previous simulation.

```r
set.seed(100)
test = simdata(10000, 5000, 2, 1, 1)
test.ash = ash(test$betahat, test$betahatsd, method = "fdr")
```

```
## [1] "normal likelihood"
```

```r
z = test$betahat/test$betahatsd
p = pchisq(z^2, df = 1, lower.tail = FALSE)
plot(qval.nonp(p, get_pi0(test.ash)), test.ash$qvalue, xlab = "nonparametric q-values (ash estimate of pi0)", 
    ylab = "ash q-values")
abline(a = 0, b = 1, col = 2, lwd = 2)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


Although these simulations are by no-means comprehensive, they show that, for these the primary
difference between ash $q$ values and qvalue $q$ values is
the different estimates of $\pi_0$.

