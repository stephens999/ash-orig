
```r
setwd("~/Documents/git/ash/paper/Rcode")
require(ashr)
```

```
## Loading required package: ashr
## Loading required package: truncnorm
```

```r
require(qvalue)
```

```
## Loading required package: qvalue
```

```r

altcol = "cyan"  #colors to use
nullcol = "blue"
nc = 40  #number of bins in histograms
ncz = 100  # number of bins in z score histograms

hh = read.table("../../data/nejm_brca_release_edit.csv", sep = ",", skip = 3)
subset = apply(hh, 1, max) < 20
hh = hh[subset, ]

labs = read.table("../../data/nejm_brca_release_edit.csv", sep = ",", skip = 1, 
    nrows = 1)
labs = 1 * (labs == "BRCA1") + 2 * (labs == "BRCA2")

hh.betahat = apply(hh[, labs == 1], 1, mean) - apply(hh[, labs == 2], 1, mean)
n1 = sum(labs == 1)
n2 = sum(labs == 2)
hh.sebetahat = sqrt(apply(hh[, labs == 1], 1, var)/n1 + apply(hh[, labs == 2], 
    1, var)/n2)
hh.zscore = hh.betahat/hh.sebetahat
hh.pval = pchisq(hh.zscore^2, df = 1, lower.tail = F)
pdf("figures/FDR_illus.pdf")
nc = 40  #number of bins in histogram
hh.hist = hist(hh.pval, freq = FALSE, xlab = "p value", main = "Distribution of p values for Hedenfalk et al data", 
    nclass = nc, col = altcol)

hh.q = qvalue(hh.pval)
abline(h = hh.q$pi0, col = nullcol, lwd = 2)

hh.hist$density = rep(hh.q$pi0, length(hh.hist$density))
# hh.hist$counts=rep(hh.q$pi0*length(hh.pval)/nc,length(hh.hist$counts))
plot(hh.hist, add = TRUE, col = nullcol, freq = FALSE)

abline(v = 0.1, lwd = 2, col = 2)

text(0.05, 1.1, labels = "A", col = 2, cex = 1.7)
text(0.05, 0.4, labels = "B", col = 2, cex = 1.7)
text(0.4, 3, labels = paste0("FDR = B/(A+B) =  ", round(hh.q$pi0 * 0.1 * length(hh.pval)/sum(hh.pval < 
    0.1), 2)), cex = 1.7)
dev.off()
```

```
## pdf 
##   2
```


Now I wanted to investigate what distribution for z values under the alternative the standard approach implicitly corresponds to.
We use three different approaches.
The first, fdrtool, we use to directly model the p values. fdrtool assumes
that the density of p values under the alternative is monotonically decreasing away from 0, and fits a non-parametric maximum likelihood distribution under this asssumption.

The other two approaches, locfdr and mixfdr, we use to directly model the z scores.

In each case, the distribution of z scores under the alternative is implied by a) the overall distribution of z scores, f(z); and b) the localfdr for each z score. Specifically $f_1(z)=(1-lfdr(z)) * f(z)$.

Note that mixfdr does not make the zero assumption, and actually models
the z score errors, so its estimated distribution for f1 is at least theoretically possible. However, we still view thie bimodal distribution as generally implausible.

Our conclusion is that making the zero assumption leads to distributions of the z scores under the
alternative that are generally implausible.

```r
# plot a histogram of z scores, highlighting the alternative distribution of
# z scores that is implied by localfdr values lfdr.
nullalthist = function(z, lfdr, ...) {
    h = hist(z, freq = FALSE, col = nullcol, nclass = ncz, ...)
    avlfdr = unlist(lapply(split(lfdr, cut(z, h$breaks), drop = FALSE), mean))
    h$density = (1 - avlfdr) * h$density
    plot(h, add = TRUE, col = altcol, freq = FALSE)
}

pdf("figures/zhist.pdf")
par(mfcol = c(3, 1))
require(fdrtool)
```

```
## Loading required package: fdrtool
```

```r
hh.fdrtool = fdrtool(hh.pval, statistic = "pvalue", plot = FALSE)
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```r
nullalthist(hh.zscore, hh.fdrtool$lfdr, main = "fdrtool")

require(locfdr)
```

```
## Loading required package: locfdr
## Loading required package: splines
```

```r
hh.locfdr = locfdr(hh.zscore, nulltype = 0, plot = 0)
nullalthist(hh.zscore, hh.locfdr$fdr, main = "locfdr")

require(mixfdr)
```

```
## Loading required package: mixfdr
```

```r
hh.mixfdr = mixFdr(hh.zscore, noiseSD = 1, theonull = TRUE, plot = FALSE)
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
## pi =		0.8016	0.0459	0.1526	
## 
## mu = 		 0.000	-2.627	 2.390	
## 
## sigma = 	1.000	1.000	1.097	
## 
## noiseSD = 	1	
```

```r
nullalthist(hh.zscore, hh.mixfdr$fdr, main = "mixfdr")
dev.off()
```

```
## pdf 
##   2
```

```r

par(mfcol = c(1, 1))
pdf("figures/zhist_fdrtool.pdf")
require(fdrtool)
hh.fdrtool = fdrtool(hh.pval, statistic = "pvalue", plot = FALSE)
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```r
nullalthist(hh.zscore, hh.fdrtool$lfdr, main = "fdrtool")
dev.off()
```

```
## pdf 
##   2
```

```r

pdf("figures/zhist_locfdr.pdf")
require(locfdr)
hh.locfdr = locfdr(hh.zscore, nulltype = 0, plot = 0)
nullalthist(hh.zscore, hh.locfdr$fdr, main = "locfdr")
dev.off()
```

```
## pdf 
##   2
```

```r

pdf("figures/zhist_mixfdr.pdf")
require(mixfdr)
hh.mixfdr = mixFdr(hh.zscore, noiseSD = 1, theonull = TRUE, plot = FALSE)
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
## pi =		0.8016	0.0459	0.1526	
## 
## mu = 		 0.000	-2.627	 2.390	
## 
## sigma = 	1.000	1.000	1.097	
## 
## noiseSD = 	1	
```

```r
nullalthist(hh.zscore, hh.mixfdr$fdr, main = "mixfdr")
dev.off()
```

```
## pdf 
##   2
```



```r
hh.ashz = ash(hh.zscore, 1, usePointMass = TRUE, VB = TRUE, prior = "nullbiased")

# now try a different prior to check robustness
mixsd = c(0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12, 10.24)
prior = rep(1, length(mixsd))
prior[1] = 100
hh.ashz2 = ash(hh.zscore, 1, mixsd = mixsd, prior = prior, VB = TRUE)
```

```
## Warning: the condition has length > 1 and only the first element will be used
## Warning: the condition has length > 1 and only the first element will be used
```

```r

pdf("figures/zhist_ash.pdf")
nullalthist(hh.zscore, hh.ashz$ZeroProb, main = "ash")
dev.off()
```

```
## pdf 
##   2
```


Compare pi0 estimates

```r
hh.fdrtool$param[3]
```

```
## [1] 0.6414
```

```r
hh.locfdr$fp0[1, 3]
```

```
## [1] 0.7432
```

```r
hh.mixfdr$pi[1]
```

```
## [1] 0.8016
```

```r
hh.ashz$fitted.g$pi[1]
```

```
## [1] 0.4658
```



```r
sum(hh.fdrtool$lfdr < 0.05)
```

```
## [1] 154
```

```r
sum(hh.locfdr$fdr < 0.05)
```

```
## [1] 171
```

```r
sum(hh.mixfdr$fdr < 0.05)
```

```
## [1] 162
```

```r
sum(hh.ashz$ZeroProb < 0.05)
```

```
## [1] 197
```






