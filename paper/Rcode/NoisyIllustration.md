The following R code simulates some data where the effects are normally distributed, and the first 500 observations have good precision, and the next 500 have poor precision.


```r
# install q value package source('http://bioconductor.org/biocLite.R')
# biocLite('qvalue')
setwd("~/Documents/git/ash/Rcode/")
library("qvalue")
library("lattice")  #library for some of the plots

# set up some data with mixture of values of s
set.seed(100)
s.good = 0.5
s.poor = 10
J.good = 500
J.poor = 500
J = J.good + J.poor
beta = c(rnorm(J, 0, 1))
s = c(rep(s.good, J.good), rep(s.poor, J.poor))
betahat = beta + rnorm(J, 0, s)
# compute the usual zscore and corresponding p value
zscore = betahat/s
pval = pchisq(zscore^2, df = 1, lower.tail = F)
```

As expected the $p$ values from the poor observations are approximately uniform, whereas those from the good observations  are enriched for small $p$ values:


```r
p = data.frame(pval = pval, type = c(rep("GOOD", 500), rep("POOR", 500)))
histogram(~pval | type, p, layout = c(1, 2), xlim = c(0, 1), breaks = seq(from = 0, 
    to = 1, length = 20), type = "count")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


Now what happens if we apply FDR methods to the all the data is
that the uniform $p$ values from the poor observations
add noise relative to looking at the good observations only. This impacts both the estimate of $pi_0$ from qvalue, and the number of findings that are significant at a given FDR.

```r
qq.all = qvalue(p$pval)
qq.good = qvalue(p$pval[p$type == "GOOD"])
print(c(qq.good$pi0, qq.all$pi0))
```

```
## [1] 0.3368 0.6932
```

```r
print(c(sum(qq.good$qvalues < 0.05), sum(qq.all$qvalues < 0.05)))
```

```
## [1] 193 114
```


In contrast, if you use adaptive shrinkage, it makes no difference whether you use all the observations or the just the good ones -
the noise is ignored (as it should be!)


```r
library("ashr")
```

```
## Loading required package: truncnorm
```

```r
beta.ash.all = ash(betahat, s)
beta.ash.good = ash(betahat[p$type == "GOOD"], s[p$type == "GOOD"])
plot(beta.ash.good$PosteriorMean, beta.ash.all$PosteriorMean[p$type == "GOOD"])
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) 

```r
plot(beta.ash.good$PositiveProb, beta.ash.all$PositiveProb[p$type == "GOOD"])
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) 


