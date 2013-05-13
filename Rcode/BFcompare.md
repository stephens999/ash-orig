try comparing BF for normal vs servin and stephens bf


```r
lABF = function(k, T) {
    return(((k/2) * T^2))
}
lBF = function(k, T, n) {
    return((-n/2) * log(1 - (k/n) * T^2))
}

n = seq(1, 100, length = 100)
T = 2
k = 0.4/T^2

plot(n, lBF(k, T, n) - lABF(k, T), type = "l")
abline(h = 0, col = 2)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 

Note that lBF>lABF always. So lABF is "biased" towards the null?
Is that right? Seems wrong, as the normal assumption is somehow
less conservative than the t?

Of course x<log(1+x)

```r
x = seq(0, 1, length = 100)
plot(x, log(1 + x))
abline(a = 0, b = 1)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


Now transform T to Z be

```r
effective.Zscore = function(T, nu) {
    return(qnorm(pt(T, nu)))
}
Z = effective.Zscore(T, n)
plot(n, lBF(k, T, n) - lABF(k, Z), type = "l")
abline(h = 0, col = 2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


