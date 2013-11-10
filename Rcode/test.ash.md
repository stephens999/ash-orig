
```r
source("ash.R")
set.seed(100)
# Test VBEM
abf = rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 1, 0), c(0, 
    0, 1, 0))
eps = 1e-10
abf[abf == 0] = eps  #replace 0 with small number
print(all.equal(VBEM(abf, c(1, 1, 1, 1))$post, c(2, 2, 4, 1)))
```

```
## [1] TRUE
```

```r
print(all.equal(VBEM(abf, c(1, 2, 1, 1))$post, c(2, 3, 4, 1)))
```

```
## [1] TRUE
```

```r

```






```r

source("ash.R")
set.seed(100)


# simulate n beta-hat values, nnull under the null with altmean and altsd
# being the mean and sd of beta under the alternative
simdata = function(n, nnull, altmean, altsd, betasd) {
    null = c(rep(1, nnull), rep(0, n - nnull))
    beta = c(rep(0, nnull), rnorm(n - nnull, altmean, altsd))
    betahat = rnorm(n, beta, betasd)
    return(list(null = null, beta = beta, betahat = betahat, betasd = betasd))
}
```



```r
ss = simdata(10000, 8000, 0, 2, 1)
beta.ash = ash(ss$betahat, ss$betasd)

system.time((post = posterior_dist(beta.ash$fitted.g, ss$betahat, rep(ss$betasd, 
    length(ss$betahat)))))
```

```
##    user  system elapsed 
##   0.128   0.008   0.136
```

```r
system.time((PositiveProb = pnormmix(0, post$pi, post$mu, post$sigma, lower.tail = FALSE)))
```

```
##    user  system elapsed 
##   0.031   0.001   0.032
```

```r

system.time((PP2 = 1 - cdf_post(beta.ash$fitted.g, 0, ss$betahat, rep(ss$betasd, 
    length(ss$betahat)))))
```

```
##    user  system elapsed 
##   0.086   0.006   0.092
```

```r
all.equal(PP2, PositiveProb)
```

```
## [1] TRUE
```




```r
system.time((beta.ash = ash(ss$betahat, ss$betasd)))
```

```
##    user  system elapsed 
##   3.075   0.087   3.164
```

```r
system.time((beta.ash.auto = ash(ss$betahat, ss$betasd, auto = TRUE)))
```

```
##    user  system elapsed 
##   1.956   0.156   2.115
```

```r
system.time((beta.ash.vb.uniform = ash(ss$betahat, ss$betasd, auto = TRUE, VB = TRUE, 
    prior = "uniform")))
```

```
##    user  system elapsed 
##  34.400   1.311  35.734
```

```r
system.time((beta.ash.vb.null = ash(ss$betahat, ss$betasd, auto = TRUE, VB = TRUE, 
    prior = NULL)))
```

```
##    user  system elapsed 
##   1.150   0.026   1.176
```

```r



hist(ss$beta, prob = TRUE, breaks = seq(-7, 7, length = 20))
x = seq(-4, 4, length = 10000)
lines(x, density(beta.ash, x), col = 2)
lines(x, density(beta.ash.auto, x), col = 3)
lines(x, density(beta.ash.vb.uniform, x), col = 4)
lines(x, density(beta.ash.vb.null, x), col = 5)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) 

```r
beta.ash$fitted.f
```

```
## NULL
```

```r
beta.ash.auto$fitted.f
```

```
## NULL
```

```r
beta.ash.vb.uniform$fitted.f
```

```
## NULL
```

```r

ss = simdata(10000, 10000, 0, 2, 1)
system.time((beta.ash = ash(ss$betahat, ss$betasd)))
```

```
##    user  system elapsed 
##   2.709   0.211   2.921
```

```r
system.time((beta.ash.auto = ash(ss$betahat, ss$betasd, auto = TRUE)))
```

```
##    user  system elapsed 
##   2.058   0.296   2.354
```

```r
system.time((beta.ash.vb.uniform = ash(ss$betahat, ss$betasd, auto = TRUE, VB = TRUE, 
    prior = "uniform")))
```

```
##    user  system elapsed 
##  26.342   3.254  29.595
```

```r
system.time((beta.ash.vb.null = ash(ss$betahat, ss$betasd, auto = TRUE, VB = TRUE, 
    prior = NULL)))
```

```
##    user  system elapsed 
##   0.290   0.033   0.323
```

```r

hist(ss$beta, prob = TRUE, breaks = seq(-7, 7, length = 20))
x = seq(-4, 4, length = 10000)
lines(x, density(beta.ash, x), col = 2)
lines(x, density(beta.ash.auto, x), col = 3)
lines(x, density(beta.ash.vb.uniform, x), col = 4)
lines(x, density(beta.ash.vb.null, x), col = 5)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) 

```r
summary(beta.ash)
```

```
## $pi
##  [1] 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## 
## $mean
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## 
## $sd
##  [1] 0.00025 0.00050 0.00100 0.00200 0.00400 0.00800 0.01600 0.03200
##  [9] 0.06400 0.12800 0.25600 0.51200 1.02400 2.04800 4.09600 8.19200
## 
## attr(,"row.names")
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
## attr(,"class")
## [1] "normalmix"
## [1] -1139.703595
## [1] TRUE
```

```r
summary(beta.ash.auto)
```

```
## $pi
## [1] 1 0 0 0 0 0 0 0
## 
## $mean
## [1] 0 0 0 0 0 0 0 0
## 
## $sd
## [1] 0.0541 0.1082 0.2164 0.4328 0.8656 1.7312 3.4623 6.9247
## 
## attr(,"row.names")
## [1] 1 2 3 4 5 6 7 8
## attr(,"class")
## [1] "normalmix"
## [1] -1131.542072
## [1] TRUE
```

```r
summary(beta.ash.vb.uniform)
```

```
## $pi
## [1] 0.5527028 0.3387469 0.0852719 0.0187021 0.0036512 0.0005927 0.0001998
## [8] 0.0001325
## 
## $mean
## [1] 0 0 0 0 0 0 0 0
## 
## $sd
## [1] 0.0541 0.1082 0.2164 0.4328 0.8656 1.7312 3.4623 6.9247
## 
## attr(,"row.names")
## [1] 1 2 3 4 5 6 7 8
## attr(,"class")
## [1] "normalmix"
## [1] -11198.50085
## [1] FALSE
```

```r
summary(beta.ash.vb.null)
```

```
## $pi
## [1] 9.999e-01 1.435e-05 1.435e-05 1.435e-05 1.435e-05 1.434e-05 1.433e-05
## [8] 1.431e-05
## 
## $mean
## [1] 0 0 0 0 0 0 0 0
## 
## $sd
## [1] 0.0541 0.1082 0.2164 0.4328 0.8656 1.7312 3.4623 6.9247
## 
## attr(,"row.names")
## [1] 1 2 3 4 5 6 7 8
## attr(,"class")
## [1] "normalmix"
## [1] -1140.81654
## [1] TRUE
```

