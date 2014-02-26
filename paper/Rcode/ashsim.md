
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





```r
pdf("figures/simABC_egdens.pdf", width = 6.5, height = 2)
plot_examples(list(simA, simB, simC))
```

```
## Warning: Removed 15 rows containing missing values (geom_path).
## Warning: Removed 15 rows containing missing values (geom_path).
## Warning: Removed 15 rows containing missing values (geom_path).
```

```r
dev.off()
```

```
## pdf 
##   2
```




```r
pdf("figures/rmse_biplot.pdf", width = 6.5, height = 2)
plot_rmse(list(simA, simB, simC))
dev.off()
```

```
## pdf 
##   2
```

```r

# This figure not used in paper?
pdf("figures/rmse_biplot_withzerobetahat.pdf")
plot_rmse(list(simA, simB, simC), inczero = TRUE, incbetahat = TRUE)
dev.off()
```

```
## pdf 
##   2
```





```r

pdf("figures/rmse_boxplot.pdf", width = 6.5, height = 2)
plot_rmse_boxplot(list(simA, simB, simC))
dev.off()
```

```
## pdf 
##   2
```

```r

pdf("figures/rmse_boxplot_extended.pdf", width = 6.5, height = 2)
plot_rmse_boxplot(list(simA, simB, simC), TRUE, TRUE, TRUE)
dev.off()
```

```
## pdf 
##   2
```

```r

```



```r

pdf("figures/loglik_boxplot.pdf", width = 6.5, height = 2)
plot_loglik_boxplot(list(simA, simB, simC))
dev.off()
```

```
## pdf 
##   2
```

```r

```





```r

pdf("figures/rmse_loglik_boxplot.pdf", width = 6.5, height = 5)
plot_rmse_loglik_boxplot(list(simA, simB, simC))
```

```
## Error: could not find function "loglik_conv"
```

```r
dev.off()
```

```
## pdf 
##   2
```


looking to see why mixfdr gives worse rmse, even though loglik looks ok.
Perhaps it overshrinks?

```r
plot(simC$fit.ash.hu[[1]]$PosteriorMean, simC$fit.mixfdr[[1]]$effectSize)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-81.png) 

```r
plot(simC$fit.ash.fdr.n[[1]]$PosteriorMean, simC$fit.mixfdr[[1]]$effectSize)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-82.png) 


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

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-101.png) 

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
lines(x, cdf.mixfdr(temp, x), col = 3)
```

```
## Error: object 'x' not found
```

```r
lines(cdf.ash(simC$fit.ash.n[[1]], x), col = 2)
```

```
## Error: object 'x' not found
```

```r
lines(cdf.ash(simC$fit.ash.hu[[1]], x), col = 2, lty = 2)
```

```
## Error: object 'x' not found
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-102.png) 





```r
plot_LR = function(sims) {
    hist(unlist(lapply(sims$fit.ash.u, get_loglik)) - unlist(lapply(sims$fit.ash.n, 
        get_loglik)), xlab = "loglik difference", main = "loglik differences for nullbiased prior vs mle", 
        nclass = 10)
}

pdf("figures/logLR.pdf")
plot_LR(simA)
plot_LR(simB)
dev.off()
```

```
## pdf 
##   2
```


## Unused figures?


```r

pdf("figures/simABC_eg_withfit.pdf", width = 6.5, height = 2)
plot_examples_withfit(list(simA, simB, simC))
```

```
## Warning: Removed 44 rows containing missing values (geom_path).
## Warning: Removed 44 rows containing missing values (geom_path).
## Warning: Removed 44 rows containing missing values (geom_path).
## Warning: Removed 44 rows containing missing values (geom_path).
## Warning: Removed 44 rows containing missing values (geom_path).
## Warning: Removed 44 rows containing missing values (geom_path).
## Warning: Removed 44 rows containing missing values (geom_path).
## Warning: Removed 44 rows containing missing values (geom_path).
## Warning: Removed 44 rows containing missing values (geom_path).
```

```r
dev.off()
```

```
## pdf 
##   2
```



```r
pdf("figures/simABC_eg_cdf_withfit.pdf", width = 6.5, height = 2)
plot_examples_cdf_withfit(list(simA, simB, simC))
```

```
## Error: no applicable method for 'mixcdf' applied to an object of class
## "normalmix"
```

```r
dev.off()
```

```
## pdf 
##   2
```


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
```

```
## Error: no applicable method for 'mixcdf' applied to an object of class
## "normalmix"
```

```r

truth = t(cdf.ash(simA$fit.ash.true[[1]], x)$y)

res.melt = melt(res, id.vars = c("beta"), variable.name = "Method")
```

```
## Error: object 'res' not found
```

```r
res.melt$Penalization = ifelse((res.melt$Method == "mixfdr" | res.melt$Method == 
    "ash.fdr.n"), "Default", "Minimal")
```

```
## Error: object 'res.melt' not found
```

```r
res.melt$Methodtype = as.character(res.melt$Method)
```

```
## Error: object 'res.melt' not found
```

```r
res.melt$Methodtype[res.melt$Method == "mixfdr.P0"] = "mixfdr"
```

```
## Error: object 'res.melt' not found
```

```r
res.melt$Methodtype[res.melt$Method == "ash.fdr.n"] = "ash"
```

```
## Error: object 'res.melt' not found
```

```r
res.melt$Methodtype[res.melt$Method == "ash.n"] = "ash"
```

```
## Error: object 'res.melt' not found
```

```r

cbbPalette <- c("#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#E69F00", 
    "#D55E00", "#CC79A7")

pdf("figures/penaltycompare.pdf")
ggplot(res.melt, aes(x = beta)) + geom_line(aes(x = beta, y = value, color = Methodtype, 
    linetype = Penalization, group = Method), size = 1.5, alpha = 0.8) + geom_line(data = res, 
    aes(y = truth, color = "truth"), alpha = 1, size = 0.5) + scale_colour_manual(name = "Method", 
    values = cbbPalette)
```

```
## Error: object 'res.melt' not found
```

```r
dev.off()
```

```
## pdf 
##   2
```


