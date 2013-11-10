

```r
source("ash.R")
beta = runif(1000)
sebetahat = 0.2
betahat = beta + rnorm(1000, 0, sebetahat)

b = seq(-5, 5, length = 100)
dd = dconvolve.uninorm.matrix(betahat, sebetahat, b)
temp = VBEM(dd, 0.1)
temp2 = b
for (i in 1:length(b)) {
    temp2[i] = dunimix(b[i], normalize(temp$post), 0, b)
}
plot(b, temp2, type = "l")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 


