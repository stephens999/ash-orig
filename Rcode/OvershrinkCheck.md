
```r
setwd("~/Documents/git/ash/Rcode/")
set.seed(32327)
## load ash functions
source("../Rcode/ash.R")
```



```r
betahat = rnorm(1000, 0, 1)
sd = rep(1, 1000)
betahat[1] = 4
hist(betahat)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-21.png) 

```r
betahat.ash = ash(betahat, sd)
betahat.ash$fitted.f
```

```
## $pi
##  [1] 0.6697469 0.0284857 0.0284857 0.0284856 0.0284855 0.0284850 0.0284831
##  [8] 0.0284753 0.0284442 0.0283205 0.0278368 0.0260186 0.0192778 0.0009694
## [15] 0.0000000 0.0000000
## 
## $sigma
##  [1] 0.00025 0.00050 0.00100 0.00200 0.00400 0.00800 0.01600 0.03200
##  [9] 0.06400 0.12800 0.25600 0.51200 1.02400 2.04800 4.09600 8.19200
## 
## $mu
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

```r
plot(betahat.ash$PosteriorMean, betahat)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-22.png) 

```r

betahat.ash2 = ash(betahat, sd, prior = "uniform")
betahat.ash2$fitted.f
```

```
## $pi
##  [1] 5.800e-01 3.625e-02 3.625e-02 3.625e-02 3.625e-02 3.625e-02 3.624e-02
##  [8] 3.623e-02 3.618e-02 3.598e-02 3.522e-02 3.244e-02 2.307e-02 3.446e-03
## [15] 2.179e-06 5.839e-12
## 
## $sigma
##  [1] 0.00025 0.00050 0.00100 0.00200 0.00400 0.00800 0.01600 0.03200
##  [9] 0.06400 0.12800 0.25600 0.51200 1.02400 2.04800 4.09600 8.19200
## 
## $mu
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

```r
plot(betahat.ash2$PosteriorMean, betahat)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-23.png) 

```r


betahat.ash3 = ash(betahat, sd, prior = 2, nullcheck = FALSE)
betahat.ash3$fitted.f
```

```
## $pi
##  [1] 0.462859 0.045543 0.045543 0.045543 0.045542 0.045542 0.045538
##  [8] 0.045524 0.045468 0.045245 0.044373 0.041092 0.029505 0.008571
## [15] 0.002592 0.001520
## 
## $sigma
##  [1] 0.00025 0.00050 0.00100 0.00200 0.00400 0.00800 0.01600 0.03200
##  [9] 0.06400 0.12800 0.25600 0.51200 1.02400 2.04800 4.09600 8.19200
## 
## $mu
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

```r
plot(betahat.ash3$PosteriorMean, betahat)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-24.png) 


Now try more extreme example


```r
betahat[1] = 8
betahat.ash3 = ash(betahat, sd, prior = 2, nullcheck = FALSE)
betahat.ash3$fitted.f
```

```
## $pi
##  [1] 0.451032 0.046982 0.046982 0.046982 0.046981 0.046980 0.046976
##  [8] 0.046961 0.046897 0.046645 0.045654 0.041838 0.027392 0.006346
## [15] 0.003131 0.002221
## 
## $sigma
##  [1] 0.00025 0.00050 0.00100 0.00200 0.00400 0.00800 0.01600 0.03200
##  [9] 0.06400 0.12800 0.25600 0.51200 1.02400 2.04800 4.09600 8.19200
## 
## $mu
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

```r
plot(betahat.ash3$PosteriorMean, betahat)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

