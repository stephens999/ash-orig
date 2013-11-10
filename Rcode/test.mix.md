
```r
############## HELPER FUNCTIONS ################################################# OUTPUT:
############## k by n matrix of normal densities
matdnorm = function(x, mu, sigma, log = FALSE) {
    k = length(mu)
    n = length(x)
    d = matrix(rep(x, rep(k, n)), nrow = k)
    return(matrix(dnorm(d, mu, sigma, log), nrow = k))
}
```



################################# TESTING #####################################

```r
setwd("~/Documents/git/ash/")
source("Rcode/mix.R")
```

```
## Loading required package: truncnorm
```

```r
source("Rcode/ash.R")
source("Rcode/ash.oldfuncs.R")
library("testthat")
temp = normalmix(c(0.5, 0.5), c(0, 0), c(1, 2))

expect_that(comp_mean(temp), equals(c(0, 0)))
expect_that(mixmean(temp), equals(0))
expect_that(mixsd(temp), equals(sqrt(2.5)))

comp_postmean(temp, 0, 1)
```

```
##      [,1]
## [1,]    0
## [2,]    0
```

```r
comp_postsd(temp, 0, 1)
```

```
##        [,1]
## [1,] 0.7071
## [2,] 0.8944
```

```r
postsd(temp, 0, 1)
```

```
## [1] 0.785
```

```r

set.seed(100)
beta = rnorm(100)
betahatse = abs(rnorm(100))
betahat = rnorm(100, beta, betahatse)

pd = posterior_dist(temp, betahat, betahatse)
mv = normmix.mv(pd)
expect_that(mv$mean, equals(postmean(temp, betahat, betahatse)))
expect_that(sqrt(mv$var), equals(postsd(temp, betahat, betahatse)))

x = seq(-5, 5, length = 100)
plot(x, dens(temp, x), type = "l")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-21.png) 

```r
all.equal(compdens(temp, x), matdnorm(x, temp$mean, temp$sd))
```

```
## [1] TRUE
```

```r

plot(x, mixcdf(temp, x), type = "l")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-22.png) 

```r

plot(x, compdens_conv(temp, x, 0.01)[1, ], type = "l")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-23.png) 

```r
plot(x, compdens_conv(temp, x, 0.1)[2, ], type = "l")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-24.png) 

```r
plot(x, dens_conv(temp, x, 0.2), type = "l")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-25.png) 

```r
plot(x, dens_conv(temp, x, 0.001), type = "l")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-26.png) 

```r
plot(x, dens_conv(temp, x, 10), type = "l")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-27.png) 





```r
temp2 = unimix(c(0.5, 0.5), c(-3, 3), c(-1, 4))
plot_dens(temp2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r

plot(x, compdens_conv(temp2, x, 0.01)[1, ], type = "l")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 

```r
plot(x, compdens_conv(temp2, x, 0.1)[2, ], type = "l")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-33.png) 

```r
plot(x, dens_conv(temp2, x, 0.2), type = "l")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-34.png) 

```r


compcdf_post(temp2, 0, c(1, 2, 3), c(10, 10, 10))
```

```
##      [,1] [,2] [,3]
## [1,]    1    1    1
## [2,]    0    0    0
```

```r
compcdf_post(temp2, -2, c(-2, 3), c(0.1, 1))
```

```
##      [,1]    [,2]
## [1,]  0.5 0.00902
## [2,]  0.0 0.00000
```

```r

plot_post_cdf(temp2, betahat = c(-2), sebetahat = 10)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-35.png) 

```r
plot_post_cdf(temp2, betahat = c(-2), sebetahat = 1, col = 2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-36.png) 



```r
etruncnorm(10, 20, 0, 1)
```

```
## [1] -Inf
```

```r
my_etruncnorm(10, 20, 0, 1)
```

```
## [1] 10.1
```

```r

etruncnorm(-20, -10, 0, 1)
```

```
## [1] -10.1
```

```r
my_etruncnorm(-20, -10, 0, 1)
```

```
## [1] -10.1
```

```r

etruncnorm(-20, -10, 5, 1)
```

```
## [1] -10.07
```

```r
my_etruncnorm(-20, -10, 5, 1)
```

```
## [1] -10.07
```

```r

etruncnorm(-30, -20, 5, 2)
```

```
## [1] -20.16
```

```r
my_etruncnorm(-30, -20, 5, 2)
```

```
## [1] -20.16
```

```r

etruncnorm(-Inf, -38, 0, 1)
```

```
## [1] -Inf
```

```r
my_etruncnorm(-Inf, -38, 0, 1)
```

```
## [1] -38
```



```r
comp_postmean(temp2, c(1, 2, 3), c(10, 10, 10))
```

```
##        [,1]   [,2]   [,3]
## [1,] -1.990 -1.987 -1.983
## [2,]  3.498  3.499  3.500
```

```r
comp_postmean(temp2, c(1, 2, 3), c(0.1, 0.1, 0.1))
```

```
##        [,1]   [,2]  [,3]
## [1,] -1.005 -1.003 -1.00
## [2,]  3.005  3.010  3.08
```


