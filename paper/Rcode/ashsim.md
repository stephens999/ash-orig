
```r
opts_knit$set(progress = TRUE, verbose = TRUE, root.dir = "~/Documents/git/ash/paper/Rcode")
require(ashr)
```

```
## Loading required package: ashr
## Loading required package: Rcpp
## Loading required package: truncnorm
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
load("ashsim_out.RData")
source("plot_examples.R")
source("plot_rmse.R")
```

```
## Loading required package: reshape2
```




looking to see why mixfdr gives worse rmse, even though loglik looks ok.
Perhaps it overshrinks?

```r
plot(simC$fit.ash.hu[[1]]$PosteriorMean, simC$fit.mixfdr[[1]]$effectSize)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r
plot(simC$fit.ash.fdr.n[[1]]$PosteriorMean, simC$fit.mixfdr[[1]]$effectSize)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 


problem seems to be that it overestimates pi0

```r
simC$fit.mixfdr[[1]]$pi
```

```
## [1] 0.71661 0.08298 0.20041
```

[1] 0.71660654 0.08298064 0.20041282



```r
temp = mixFdr(simC$betahat[[1]]/simC$betahatsd[[1]], noiseSD = 1, theonull = TRUE, 
    calibrate = TRUE)
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Calibrating penalization (slow) 
## To avoid this step next time, note the calibrated penalization parameter and supply it to the function. 
## Calibration: Generating random data...done
## Calibration: Testing 20 different penalization, each B = 25 times
## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20  done 
## Calibrated penalization is 100 
## Fitting final model
```

```
## Warning: Null proportion pi0 is small. Consider increasing penalization and/or using an empirical null.
## Warning: Using an empirical null with a fitted noiseSD gives a substantially different model. Consider rerunning with theonull = FALSE and noiseSD = NA.
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-51.png) 

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.6418	0.2483	0.1099	
## 
## mu = 		 0.000	-2.822	 2.306	
## 
## sigma = 	1.000	1.602	1.000	
## 
## noiseSD = 	1	
```

```r
temp2 = mixFdr(simC$betahat[[1]]/simC$betahatsd[[1]], noiseSD = 1, theonull = TRUE, 
    J = 10, plot = FALSE)
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model
```

```
## Warning: Null proportion pi0 is small. Consider increasing penalization and/or using an empirical null.
## Warning: Using an empirical null with a fitted noiseSD gives a substantially different model. Consider rerunning with theonull = FALSE and noiseSD = NA.
```

```
## 
## Fitted Model: J = 10 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	
## 
## pi =		0.7065	0.0139	0.0178	0.0123	0.0322	0.0164	0.0076	0.1549	0.0169	0.0214	
## 
## mu = 		 0.000	 2.408	 2.412	 2.407	-2.364	 2.410	 2.400	-2.904	-6.065	 2.414	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.152	1.248	1.000	
## 
## noiseSD = 	1	
```

```r
plot(ecdf(simC$beta[[1]]), xlim = c(-6, 6))
x = seq(-6, 6, length = 100)
lines(x, cdf.mixfdr(temp, x), col = 3)
lines(cdf.ash(simC$fit.ash.n[[1]], x), col = 2)
lines(cdf.ash(simC$fit.ash.hu[[1]], x), col = 2, lty = 2)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-52.png) 





Experiment with mixFdr to see how penalization effects fit. 

```r
betahat = simA$betahat[[1]]
betahatsd = simA$betahatsd[[1]]
set.seed(111)  #to ensure reproducibility
fit.mixfdr.P0 = mixFdr(betahat/betahatsd, noiseSD = 1, theonull = TRUE, P = 0, 
    plot = FALSE)
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model
```

```
## Warning: Null proportion pi0 is small. Consider increasing penalization and/or using an empirical null.
## Warning: Using an empirical null with a fitted noiseSD gives a substantially different model. Consider rerunning with theonull = FALSE and noiseSD = NA.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.3430	0.6113	0.0457	
## 
## mu = 		 0.0000	 0.2683	-2.2509	
## 
## sigma = 	1.000	1.942	1.614	
## 
## noiseSD = 	1	
```

```r
cdf.mixfdr = function(a, x) {
    mixcdf(normalmix(a$pi, a$mu, a$sigma - 1), x)
}
x = seq(-4, 4, length = 100)

res = data.frame(beta = x, mixfdr.P0 = t(cdf.mixfdr(fit.mixfdr.P0, x)), mixfdr = t(cdf.mixfdr(simA$fit.mixfdr[[1]], 
    x)), ash.n = t(cdf.ash(simA$fit.ash.n[[1]], x)$y), ash.fdr.n = t(cdf.ash(simA$fit.ash.fdr.n[[1]], 
    x)$y))

truth = t(cdf.ash(simA$fit.ash.true[[1]], x)$y)

res.melt = melt(res, id.vars = c("beta"), variable.name = "Method")
res.melt$Penalization = ifelse((res.melt$Method == "mixfdr" | res.melt$Method == 
    "ash.fdr.n"), "Default", "Minimal")
res.melt$Methodtype = as.character(res.melt$Method)
res.melt$Methodtype[res.melt$Method == "mixfdr.P0"] = "mixfdr"
res.melt$Methodtype[res.melt$Method == "ash.fdr.n"] = "ash"
res.melt$Methodtype[res.melt$Method == "ash.n"] = "ash"

cbbPalette <- c("#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#E69F00", 
    "#D55E00", "#CC79A7")

# pdf('figures/penaltycompare.pdf')
ggplot(res.melt, aes(x = beta)) + geom_line(aes(x = beta, y = value, color = Methodtype, 
    linetype = Penalization, group = Method), size = 1.5, alpha = 0.8) + geom_line(data = res, 
    aes(y = truth, color = "truth"), alpha = 1, size = 0.5) + scale_colour_manual(name = "Method", 
    values = cbbPalette)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

```r
# dev.off()
```


