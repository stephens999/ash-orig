
```r
opts_knit$set(progress = TRUE, verbose = TRUE, root.dir = "~/Documents/git/ash/paper/Rcode")
require(ashr)
```

```
## Loading required package: ashr
## Loading required package: Rcpp
```

```r
require(qvalue)
```

```
## Loading required package: qvalue
```

```r
require(fdrtool)
```

```
## Loading required package: fdrtool
```

```r
require(mixfdr)
```

```
## Loading required package: mixfdr
```

```r
require(locfdr)
```

```
## Loading required package: locfdr
## Loading required package: splines
```

```r
require(ggplot2)
```

```
## Loading required package: ggplot2
## 
## Attaching package: 'ggplot2'
## 
## The following object is masked from 'package:qvalue':
## 
##     qplot
```



```r
load("sim1_out.RData")
```


#look into methods for estimating pi0 conservatively



```r
# par(mfcol=c(3,3))
plot_ecdf = function(sims) {
    for (i in 1:length(sims$beta)) {
        plot(ecdf(sims$beta[[i]]), xlim = c(-6, 6), main = paste0("iteration ", 
            i))
        x = seq(-6, 6, length = 1000)
        lines(cdf.ash(sims$betahat.ash.n[[i]], x), col = 2, lwd = 2)
        lines(cdf.ash(sims$betahat.ash.u[[i]], x), col = 3, lwd = 2)
        lines(cdf.ash(sims$betahat.ash.true[[i]], x), col = 4, lwd = 2)
    }
}
plot_ecdf(simres1)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-33.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-34.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-35.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-36.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-37.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-38.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-39.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-310.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-311.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-312.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-313.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-314.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-315.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-316.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-317.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-318.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-319.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-320.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-321.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-322.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-323.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-324.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-325.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-326.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-327.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-328.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-329.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-330.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-331.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-332.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-333.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-334.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-335.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-336.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-337.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-338.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-339.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-340.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-341.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-342.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-343.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-344.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-345.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-346.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-347.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-348.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-349.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-350.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-351.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-352.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-353.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-354.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-355.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-356.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-357.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-358.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-359.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-360.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-361.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-362.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-363.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-364.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-365.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-366.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-367.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-368.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-369.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-370.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-371.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-372.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-373.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-374.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-375.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-376.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-377.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-378.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-379.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-380.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-381.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-382.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-383.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-384.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-385.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-386.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-387.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-388.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-389.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-390.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-391.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-392.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-393.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-394.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-395.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-396.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-397.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-398.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-399.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3100.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3101.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3102.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3103.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3104.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3105.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3106.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3107.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3108.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3109.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3110.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3111.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3112.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3113.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3114.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3115.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3116.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3117.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3118.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3119.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3120.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3121.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3122.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3123.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3124.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3125.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3126.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3127.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3128.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3129.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3130.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3131.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3132.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3133.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3134.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3135.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3136.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3137.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3138.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3139.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3140.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3141.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3142.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3143.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3144.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3145.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3146.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3147.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3148.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3149.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3150.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3151.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3152.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3153.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3154.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3155.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3156.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3157.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3158.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3159.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3160.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3161.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3162.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3163.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3164.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3165.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3166.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3167.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3168.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3169.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3170.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3171.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3172.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3173.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3174.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3175.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3176.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3177.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3178.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3179.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3180.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3181.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3182.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3183.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3184.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3185.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3186.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3187.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3188.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3189.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3190.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3191.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3192.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3193.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3194.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3195.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3196.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3197.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3198.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3199.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3200.png) 

```r
plot_ecdf(simres2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3201.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3202.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3203.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3204.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3205.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3206.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3207.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3208.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3209.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3210.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3211.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3212.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3213.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3214.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3215.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3216.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3217.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3218.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3219.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3220.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3221.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3222.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3223.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3224.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3225.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3226.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3227.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3228.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3229.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3230.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3231.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3232.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3233.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3234.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3235.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3236.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3237.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3238.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3239.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3240.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3241.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3242.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3243.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3244.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3245.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3246.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3247.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3248.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3249.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3250.png) 


Figure to show that estimated betahats are not so different

```r
plot(betahat.ash.u[[1]]$PosteriorMean, betahat.ash.n[[1]]$PosteriorMean)
```

```
## Error: object 'betahat.ash.u' not found
```

```r
abline(a = 0, b = 1, lwd = 2, col = 2)
```

```
## Error: plot.new has not been called yet
```



QUestion: is the null-biased prior maybe a little too conservative?
Answer: log likelihoods don't suggest they are



```r
# hh.ashtrue = hh.ashz hh.ashtrue$fitted.g$pi =
# c(2/3,1/15,1/15,1/15,1/15,1/15) hh.ashtrue$fitted.g$mean = c(0,0,0,0,0,0)
# hh.ashtrue$fitted.g$sd = sqrt(c(0,1,0.2,0.4,0.8,3))

# loglik(hh.ashtrue,betahat,sebetahat) loglik(hh.ashz,betahat,sebetahat)
```


