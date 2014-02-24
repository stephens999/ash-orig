#look into methods for estimating pi0 conservatively

Do some simulations


```r
opts_knit$set(progress = TRUE, verbose = TRUE, root.dir = "~/Documents/git/ash/paper/Rcode")
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
# bsd gives standard deviation of beta pi0 is simulated to be uniform on
# [minpi0,1]
basicsim = function(mixsd, mixpi_alt, bsd = 1, minpi0 = 0, seedval = 100, nsamp = 1000, 
    niter = 50) {
    set.seed(seedval)
    beta = list()
    betahatsd = list()
    betahat = list()
    zscore = list()
    pval = list()
    betahat.ash.n = list()
    betahat.ash.u = list()
    betahat.ash.npm = list()
    betahat.ash.true = list()
    betahat.qval = list()
    betahat.fdrtool = list()
    betahat.locfdr = list()
    betahat.mixfdr = list()
    pi0 = rep(0, niter)
    for (i in 1:niter) {
        pi0[i] = runif(1, minpi0, 1)
        mixpi = c(pi0[i], (1 - pi0[i]) * mixpi_alt)
        sd = sample(mixsd, nsamp, prob = mixpi, replace = TRUE)
        beta[[i]] = rnorm(nsamp, 0, sd)
        betahatsd[[i]] = bsd
        betahat[[i]] = beta[[i]] + rnorm(nsamp, 0, betahatsd[[i]])
        zscore[[i]] = betahat[[i]]/betahatsd[[i]]
        pval[[i]] = pchisq(zscore[[i]]^2, df = 1, lower.tail = F)
        betahat.ash.n[[i]] = ash(betahat[[i]], betahatsd[[i]], pointmass = TRUE, 
            prior = "nullbiased", gridmult = 2)
        betahat.ash.u[[i]] = ash(betahat[[i]], betahatsd[[i]], pointmass = TRUE, 
            prior = "uniform", gridmult = 2)
        betahat.ash.npm[[i]] = ash(betahat[[i]], betahatsd[[i]], pointmass = FALSE, 
            prior = "uniform", gridmult = 2)
        betahat.ash.true[[i]] = ash(betahat[[i]], betahatsd[[i]], g = normalmix(mixpi, 
            rep(0, length(mixpi)), mixsd))
        
        betahat.qval[[i]] = qvalue(pval[[i]])
        betahat.fdrtool[[i]] = fdrtool(pval[[i]], statistic = "pvalue", plot = FALSE)
        betahat.locfdr[[i]] = locfdr(zscore[[i]], nulltype = 0, plot = 0)
        betahat.mixfdr[[i]] = mixFdr(zscore[[i]], noiseSD = 1, theonull = TRUE, 
            plot = FALSE)
    }
    return(list(beta = beta, betahatsd = betahatsd, betahat = betahat, zscore = zscore, 
        pval = pval, betahat.ash.n = betahat.ash.n, betahat.ash.u = betahat.ash.u, 
        betahat.ash.npm = betahat.ash.npm, betahat.ash.true = betahat.ash.true, 
        betahat.qval = betahat.qval, betahat.fdrtool = betahat.fdrtool, betahat.locfdr = betahat.locfdr, 
        betahat.mixfdr = betahat.mixfdr, pi0 = pi0))
}
```



```r
mixsd = c(0, 0.25, 0.5, 1, 2)
mixpi_alt = c(0.4, 0.2, 0.2, 0.2)  #mixture proportions under the alternative

simres1 = basicsim(mixsd, mixpi_alt, niter = 200, nsamp = 1000)
```

```
## Loading required package: fdrtool
## Loading required package: mixfdr
## Loading required package: locfdr
## Loading required package: splines
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9415	0.0215	0.0371	
## 
## mu = 		 0.000	-2.903	 2.983	
## 
## sigma = 	1.000	1.348	1.093	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9538	0.0200	0.0262	
## 
## mu = 		 0.000	 3.065	-3.279	
## 
## sigma = 	1.000	1.507	1.199	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9557	0.0325	0.0117	
## 
## mu = 		 0.000	-2.394	 3.255	
## 
## sigma = 	1.000	1.000	1.481	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9611	0.0175	0.0215	
## 
## mu = 		 0.000	-3.053	 3.404	
## 
## sigma = 	1.000	1.255	1.720	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8991	0.0529	0.0481	
## 
## mu = 		 0.000	 2.724	-2.924	
## 
## sigma = 	1.000	1.124	1.091	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9790	0.0113	0.0097	
## 
## mu = 		 0.000	 2.630	-2.916	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9555	0.0222	0.0223	
## 
## mu = 		 0.000	 3.141	-2.879	
## 
## sigma = 	1.000	1.000	1.085	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9577	0.0203	0.0220	
## 
## mu = 		 0.000	-3.017	 3.480	
## 
## sigma = 	1.000	1.590	1.337	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9839	0.0085	0.0076	
## 
## mu = 		 0.000	 2.900	-3.417	
## 
## sigma = 	1.000	1.000	1.001	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9724	0.0155	0.0121	
## 
## mu = 		 0.000	-3.059	 2.982	
## 
## sigma = 	1.000	1.159	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9388	0.0323	0.0289	
## 
## mu = 		 0.000	-2.754	 3.300	
## 
## sigma = 	1.000	1.193	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9260	0.0349	0.0391	
## 
## mu = 		 0.000	-3.115	 2.802	
## 
## sigma = 	1.000	1.182	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9532	0.0130	0.0338	
## 
## mu = 		 0.000	-3.668	 3.145	
## 
## sigma = 	1.000	1.597	1.302	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9714	0.0109	0.0177	
## 
## mu = 		 0.000	-2.721	 2.741	
## 
## sigma = 	1.000	1.503	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9969	0.0019	0.0012	
## 
## mu = 		 0.000	 2.942	-1.582	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9782	0.0071	0.0147	
## 
## mu = 		 0.000	 3.043	-2.927	
## 
## sigma = 	1.000	1.063	1.026	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9939	0.0051	0.0010	
## 
## mu = 		 0.000	-2.661	 4.777	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9454	0.0263	0.0283	
## 
## mu = 		 0.000	-3.120	 2.732	
## 
## sigma = 	1.000	1.141	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9702	0.0142	0.0155	
## 
## mu = 		 0.000	 2.934	-3.007	
## 
## sigma = 	1.000	1.722	1.275	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9529	0.0269	0.0202	
## 
## mu = 		 0.000	 2.978	-3.445	
## 
## sigma = 	1.000	1.287	1.123	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9195	0.0441	0.0364	
## 
## mu = 		 0.000	 2.816	-3.097	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9351	0.0417	0.0232	
## 
## mu = 		 0.000	 2.636	-3.357	
## 
## sigma = 	1.000	1.000	1.355	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9622	0.0192	0.0186	
## 
## mu = 		 0.000	-3.162	 2.876	
## 
## sigma = 	1.000	1.246	1.081	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9397	0.0411	0.0192	
## 
## mu = 		 0.000	 2.769	-2.945	
## 
## sigma = 	1.000	1.002	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9732	0.0134	0.0134	
## 
## mu = 		 0.000	 3.312	-2.945	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9708	0.0204	0.0088	
## 
## mu = 		 0.000	-2.822	 3.111	
## 
## sigma = 	1.000	1.255	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9711	0.0200	0.0090	
## 
## mu = 		 0.000	 3.288	-3.149	
## 
## sigma = 	1.000	1.512	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9811	0.0061	0.0128	
## 
## mu = 		 0.000	 3.024	-3.304	
## 
## sigma = 	1.000	1.000	1.194	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9848	0.0074	0.0078	
## 
## mu = 		 0.000	-3.092	 2.293	
## 
## sigma = 	1.000	1.196	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9807	0.0095	0.0098	
## 
## mu = 		 0.000	-3.299	 2.603	
## 
## sigma = 	1.000	1.143	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9479	0.0248	0.0274	
## 
## mu = 		 0.000	-2.994	 3.158	
## 
## sigma = 	1.000	1.272	1.593	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9808	0.0104	0.0088	
## 
## mu = 		 0.000	-2.661	 3.815	
## 
## sigma = 	1.000	1.064	1.483	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9832	0.0121	0.0047	
## 
## mu = 		 0.000	 3.765	-4.055	
## 
## sigma = 	1.000	1.000	1.478	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9792	0.0082	0.0126	
## 
## mu = 		 0.000	-3.032	 3.372	
## 
## sigma = 	1.000	1.000	1.154	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9708	0.0195	0.0097	
## 
## mu = 		 0.000	 3.059	-3.392	
## 
## sigma = 	1.000	1.507	1.202	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Using an empirical null with a fitted noiseSD gives a
## substantially different model. Consider rerunning with theonull = FALSE
## and noiseSD = NA.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9211	0.0440	0.0350	
## 
## mu = 		 0.000	 2.768	-3.047	
## 
## sigma = 	1.000	1.355	1.384	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9380	0.0349	0.0271	
## 
## mu = 		 0.000	 2.990	-2.905	
## 
## sigma = 	1.000	1.115	1.318	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9889	0.0044	0.0068	
## 
## mu = 		 0.000	-3.240	 3.304	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9545	0.0235	0.0219	
## 
## mu = 		 0.000	-2.683	 3.458	
## 
## sigma = 	1.000	1.000	1.446	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9917	0.0016	0.0067	
## 
## mu = 		 0.000	-1.560	 2.802	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9651	0.0117	0.0232	
## 
## mu = 		 0.000	 3.371	-3.207	
## 
## sigma = 	1.000	1.117	1.455	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9861	0.0086	0.0052	
## 
## mu = 		 0.000	-2.806	 2.687	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9750	0.0176	0.0074	
## 
## mu = 		 0.000	 2.489	-3.579	
## 
## sigma = 	1.000	1.000	1.599	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9603	0.0206	0.0191	
## 
## mu = 		 0.000	 2.865	-3.114	
## 
## sigma = 	1.000	1.000	1.114	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9500	0.0248	0.0252	
## 
## mu = 		 0.000	 2.016	-2.630	
## 
## sigma = 	1.000	2.283	1.161	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9579	0.0211	0.0209	
## 
## mu = 		 0.000	 3.137	-3.176	
## 
## sigma = 	1.000	1.317	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9406	0.0243	0.0351	
## 
## mu = 		 0.000	 2.983	-2.846	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9755	0.0163	0.0081	
## 
## mu = 		 0.000	 2.673	-3.246	
## 
## sigma = 	1.000	1.000	1.147	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.944	0.034	0.022	
## 
## mu = 		 0.000	-2.592	 3.549	
## 
## sigma = 	1.000	1.247	1.474	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9262	0.0354	0.0384	
## 
## mu = 		 0.000	 2.868	-2.857	
## 
## sigma = 	1.000	1.478	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9959	0.0011	0.0030	
## 
## mu = 		 0.000	 0.978	-3.837	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9810	0.0142	0.0048	
## 
## mu = 		 0.000	-2.397	 4.423	
## 
## sigma = 	1.000	1.000	2.399	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	 TRUE	FALSE	
## 
## pi =		0.9983	0.0010	0.0007	
## 
## mu = 		 0.0000	-0.0783	-3.3539	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9681	0.0152	0.0168	
## 
## mu = 		 0.000	-3.057	 3.402	
## 
## sigma = 	1.000	1.247	1.237	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Using an empirical null with a fitted noiseSD gives a
## substantially different model. Consider rerunning with theonull = FALSE
## and noiseSD = NA.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9179	0.0353	0.0469	
## 
## mu = 		 0.000	-2.840	 2.645	
## 
## sigma = 	1.000	1.037	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9498	0.0196	0.0306	
## 
## mu = 		 0.000	 3.059	-2.780	
## 
## sigma = 	1.000	1.019	1.283	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9387	0.0271	0.0342	
## 
## mu = 		 0.000	 2.991	-3.145	
## 
## sigma = 	1.000	1.025	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9966	0.0015	0.0019	
## 
## mu = 		 0.000	-2.864	 3.398	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9725	0.0151	0.0124	
## 
## mu = 		 0.000	-2.487	 3.172	
## 
## sigma = 	1.000	1.000	1.025	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9575	0.0214	0.0210	
## 
## mu = 		 0.000	 2.895	-3.148	
## 
## sigma = 	1.000	1.487	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9713	0.0133	0.0154	
## 
## mu = 		 0.000	 2.787	-3.038	
## 
## sigma = 	1.000	1.524	1.455	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9637	0.0210	0.0153	
## 
## mu = 		 0.000	 2.752	-3.361	
## 
## sigma = 	1.000	1.000	1.237	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9384	0.0260	0.0356	
## 
## mu = 		 0.000	-2.885	 2.748	
## 
## sigma = 	1.000	1.303	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9266	0.0347	0.0387	
## 
## mu = 		 0.000	-2.714	 2.864	
## 
## sigma = 	1.000	1.315	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9948	0.0021	0.0032	
## 
## mu = 		 0.000	 4.457	-2.420	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9809	0.0095	0.0096	
## 
## mu = 		 0.000	-3.555	 3.221	
## 
## sigma = 	1.000	1.000	1.353	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.961	0.015	0.024	
## 
## mu = 		 0.000	 3.010	-3.143	
## 
## sigma = 	1.000	1.197	1.056	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9271	0.0284	0.0444	
## 
## mu = 		 0.000	-2.887	 2.463	
## 
## sigma = 	1.000	1.266	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Using an empirical null with a fitted noiseSD gives a
## substantially different model. Consider rerunning with theonull = FALSE
## and noiseSD = NA.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9051	0.0454	0.0495	
## 
## mu = 		 0.000	 2.585	-3.210	
## 
## sigma = 	1.000	1.145	1.158	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9109	0.0362	0.0529	
## 
## mu = 		 0.000	 3.284	-2.989	
## 
## sigma = 	1.000	1.006	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9556	0.0218	0.0227	
## 
## mu = 		 0.000	-2.905	 2.877	
## 
## sigma = 	1.000	1.295	1.227	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9766	0.0045	0.0189	
## 
## mu = 		 0.000	-3.541	 2.883	
## 
## sigma = 	1.000	1.322	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9880	0.0036	0.0085	
## 
## mu = 		0.000	2.887	3.189	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9892	0.0036	0.0072	
## 
## mu = 		 0.000	-3.404	 3.469	
## 
## sigma = 	1.000	1.415	1.859	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9382	0.0331	0.0287	
## 
## mu = 		 0.000	-2.973	 2.946	
## 
## sigma = 	1.000	1.393	1.442	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9967	0.0010	0.0023	
## 
## mu = 		 0.000	 6.207	-2.562	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9154	0.0579	0.0267	
## 
## mu = 		 0.000	 3.170	-3.164	
## 
## sigma = 	1.000	1.000	1.017	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9307	0.0392	0.0301	
## 
## mu = 		 0.000	 2.790	-2.554	
## 
## sigma = 	1.000	1.125	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9969	0.0018	0.0014	
## 
## mu = 		 0.000	 2.143	-1.274	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9911	0.0056	0.0032	
## 
## mu = 		0.000	2.207	2.117	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9866	0.0060	0.0074	
## 
## mu = 		 0.00	-2.90	 3.41	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9821	0.0142	0.0037	
## 
## mu = 		 0.000	 2.845	-4.000	
## 
## sigma = 	1.00	1.00	1.18	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9877	0.0068	0.0055	
## 
## mu = 		 0.000	 2.721	-2.824	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9945	0.0030	0.0025	
## 
## mu = 		0.0000	2.2668	0.8778	
## 
## sigma = 	1.000	3.572	3.453	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9750	0.0094	0.0156	
## 
## mu = 		 0.000	-3.652	 3.160	
## 
## sigma = 	1.000	1.411	1.377	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9942	0.0040	0.0018	
## 
## mu = 		0.0000	1.2871	0.4115	
## 
## sigma = 	1.000	4.283	1.385	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9823	0.0059	0.0118	
## 
## mu = 		 0.000	-3.059	 3.261	
## 
## sigma = 	1.000	1.000	1.232	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9916	0.0033	0.0051	
## 
## mu = 		 0.000	 3.370	-3.434	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9455	0.0350	0.0195	
## 
## mu = 		 0.000	 2.720	-2.691	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9755	0.0170	0.0074	
## 
## mu = 		 0.000	 3.595	-3.456	
## 
## sigma = 	1.000	1.204	1.261	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9823	0.0104	0.0072	
## 
## mu = 		 0.000	 2.939	-3.024	
## 
## sigma = 	1.000	1.000	1.127	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9748	0.0117	0.0135	
## 
## mu = 		 0.000	-3.334	 2.721	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9827	0.0070	0.0103	
## 
## mu = 		 0.000	 3.577	-3.805	
## 
## sigma = 	1.000	2.408	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9624	0.0182	0.0194	
## 
## mu = 		 0.000	-2.681	 3.093	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9945	0.0020	0.0035	
## 
## mu = 		 0.000	-5.404	 2.826	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9459	0.0202	0.0339	
## 
## mu = 		 0.000	 3.163	-2.994	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9618	0.0130	0.0251	
## 
## mu = 		 0.000	-2.778	 3.407	
## 
## sigma = 	1.000	1.088	1.253	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9745	0.0129	0.0125	
## 
## mu = 		 0.000	-2.922	 3.375	
## 
## sigma = 	1.000	1.684	1.422	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9778	0.0103	0.0120	
## 
## mu = 		 0.000	-3.447	 3.061	
## 
## sigma = 	1.000	2.021	1.425	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9603	0.0212	0.0185	
## 
## mu = 		 0.000	 2.476	-2.891	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9391	0.0272	0.0337	
## 
## mu = 		 0.000	-3.484	 2.990	
## 
## sigma = 	1.000	1.641	1.362	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9699	0.0094	0.0207	
## 
## mu = 		 0.000	 3.530	-3.279	
## 
## sigma = 	1.000	1.169	1.128	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9798	0.0100	0.0102	
## 
## mu = 		 0.000	 2.980	-3.135	
## 
## sigma = 	1.000	1.000	1.017	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9215	0.0382	0.0404	
## 
## mu = 		 0.000	 3.051	-2.865	
## 
## sigma = 	1.000	1.109	1.071	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9718	0.0062	0.0220	
## 
## mu = 		 0.000	-3.166	 2.811	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	 TRUE	FALSE	
## 
## pi =		0.9929	0.0010	0.0061	
## 
## mu = 		0.0000	0.1216	3.3177	
## 
## sigma = 	1.000	1.000	1.097	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9778	0.0074	0.0148	
## 
## mu = 		 0.000	-2.514	 3.325	
## 
## sigma = 	1.00	1.00	1.31	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9728	0.0204	0.0068	
## 
## mu = 		 0.000	-2.811	 3.151	
## 
## sigma = 	1.000	1.264	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9412	0.0233	0.0354	
## 
## mu = 		 0.000	-3.014	 3.249	
## 
## sigma = 	1.000	1.115	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9837	0.0030	0.0133	
## 
## mu = 		0.000	2.753	2.650	
## 
## sigma = 	1.000	1.107	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9678	0.0139	0.0182	
## 
## mu = 		 0.000	 2.995	-2.969	
## 
## sigma = 	1.000	1.286	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9626	0.0132	0.0242	
## 
## mu = 		 0.000	-2.996	 2.837	
## 
## sigma = 	1.000	1.423	1.021	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9451	0.0229	0.0320	
## 
## mu = 		 0.000	-3.187	 2.954	
## 
## sigma = 	1.000	1.181	1.437	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9834	0.0061	0.0106	
## 
## mu = 		 0.000	 2.996	-2.736	
## 
## sigma = 	1.000	2.003	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9425	0.0213	0.0362	
## 
## mu = 		 0.000	 3.481	-2.939	
## 
## sigma = 	1.000	1.226	1.105	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9912	0.0076	0.0011	
## 
## mu = 		 0.000	-2.996	 1.014	
## 
## sigma = 	1.000	1.199	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9919	0.0020	0.0061	
## 
## mu = 		 0.000	-5.185	 3.683	
## 
## sigma = 	1.000	1.000	1.097	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9778	0.0109	0.0113	
## 
## mu = 		 0.000	-3.268	 2.557	
## 
## sigma = 	1.000	1.214	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9132	0.0517	0.0351	
## 
## mu = 		 0.000	 2.871	-2.804	
## 
## sigma = 	1.000	1.239	1.473	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9393	0.0356	0.0252	
## 
## mu = 		 0.000	-2.878	 3.171	
## 
## sigma = 	1.000	1.259	1.208	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9376	0.0438	0.0186	
## 
## mu = 		 0.000	-2.801	 3.197	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9622	0.0178	0.0199	
## 
## mu = 		 0.000	-2.861	 3.339	
## 
## sigma = 	1.000	1.000	1.534	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9549	0.0218	0.0232	
## 
## mu = 		 0.000	-2.868	 2.775	
## 
## sigma = 	1.000	1.000	1.039	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9642	0.0162	0.0197	
## 
## mu = 		 0.000	-2.834	 2.863	
## 
## sigma = 	1.000	1.116	1.090	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9464	0.0310	0.0226	
## 
## mu = 		 0.000	 2.890	-2.669	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9874	0.0091	0.0034	
## 
## mu = 		 0.000	 3.473	-3.207	
## 
## sigma = 	1.000	1.000	1.125	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9892	0.0064	0.0044	
## 
## mu = 		 0.000	-3.773	 3.982	
## 
## sigma = 	1.000	1.000	1.831	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9332	0.0318	0.0350	
## 
## mu = 		 0.000	-3.045	 2.888	
## 
## sigma = 	1.000	1.143	1.003	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Using an empirical null with a fitted noiseSD gives a
## substantially different model. Consider rerunning with theonull = FALSE
## and noiseSD = NA.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9023	0.0500	0.0477	
## 
## mu = 		 0.000	 3.108	-2.790	
## 
## sigma = 	1.000	1.373	1.567	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9471	0.0274	0.0255	
## 
## mu = 		 0.000	-3.434	 2.729	
## 
## sigma = 	1.000	1.001	1.053	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9950	0.0019	0.0031	
## 
## mu = 		 0.000	-2.499	-2.646	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9951	0.0028	0.0021	
## 
## mu = 		 0.000	-4.137	 3.900	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9645	0.0141	0.0214	
## 
## mu = 		 0.000	-3.064	 2.670	
## 
## sigma = 	1.000	1.529	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9502	0.0211	0.0287	
## 
## mu = 		 0.000	-3.202	 2.971	
## 
## sigma = 	1.000	1.306	1.015	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9710	0.0165	0.0125	
## 
## mu = 		 0.000	-3.113	 3.008	
## 
## sigma = 	1.000	1.422	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9709	0.0106	0.0184	
## 
## mu = 		 0.000	-3.154	 2.755	
## 
## sigma = 	1.000	1.106	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9444	0.0311	0.0245	
## 
## mu = 		 0.000	 3.582	-2.892	
## 
## sigma = 	1.000	1.194	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9457	0.0325	0.0218	
## 
## mu = 		 0.000	-3.063	 2.998	
## 
## sigma = 	1.000	1.023	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9675	0.0174	0.0150	
## 
## mu = 		 0.000	-2.744	 3.837	
## 
## sigma = 	1.000	1.091	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9927	0.0024	0.0050	
## 
## mu = 		 0.000	 2.292	-2.804	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9877	0.0044	0.0079	
## 
## mu = 		 0.000	 2.904	-3.694	
## 
## sigma = 	1.00	1.04	1.00	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9174	0.0395	0.0431	
## 
## mu = 		 0.000	 2.645	-2.642	
## 
## sigma = 	1.000	1.000	1.152	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9718	0.0135	0.0148	
## 
## mu = 		 0.000	-3.483	 2.892	
## 
## sigma = 	1.000	1.153	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9511	0.0250	0.0240	
## 
## mu = 		 0.000	 2.778	-2.802	
## 
## sigma = 	1.000	1.000	1.201	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9891	0.0019	0.0090	
## 
## mu = 		 0.000	-5.083	 3.563	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9278	0.0282	0.0440	
## 
## mu = 		 0.000	-2.867	 2.964	
## 
## sigma = 	1.000	1.115	1.294	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9834	0.0077	0.0089	
## 
## mu = 		 0.000	 2.731	-3.335	
## 
## sigma = 	1.000	1.000	1.272	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9354	0.0230	0.0417	
## 
## mu = 		 0.000	-3.015	 2.530	
## 
## sigma = 	1.000	1.456	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9775	0.0138	0.0087	
## 
## mu = 		 0.000	 3.394	-3.526	
## 
## sigma = 	1.000	1.034	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9560	0.0251	0.0190	
## 
## mu = 		 0.000	 2.973	-2.960	
## 
## sigma = 	1.000	1.227	1.084	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9700	0.0138	0.0162	
## 
## mu = 		 0.000	-2.652	 2.787	
## 
## sigma = 	1.000	1.255	1.397	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9877	0.0042	0.0081	
## 
## mu = 		 0.000	-4.178	 3.056	
## 
## sigma = 	1.000	1.000	1.388	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9153	0.0528	0.0319	
## 
## mu = 		 0.000	-2.963	 3.072	
## 
## sigma = 	1.000	1.215	1.125	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9286	0.0309	0.0406	
## 
## mu = 		 0.000	-2.678	 2.735	
## 
## sigma = 	1.000	1.000	1.554	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9703	0.0142	0.0156	
## 
## mu = 		 0.000	 2.982	-2.956	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9892	0.0076	0.0032	
## 
## mu = 		 0.000	 3.047	-2.686	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9709	0.0161	0.0130	
## 
## mu = 		 0.000	-2.529	 2.488	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9383	0.0309	0.0308	
## 
## mu = 		 0.000	 3.174	-3.154	
## 
## sigma = 	1.000	1.353	1.167	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9454	0.0309	0.0236	
## 
## mu = 		 0.000	-2.780	 3.206	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9825	0.0113	0.0062	
## 
## mu = 		 0.000	-3.812	 3.558	
## 
## sigma = 	1.000	1.363	1.644	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9309	0.0348	0.0344	
## 
## mu = 		 0.000	 2.827	-3.098	
## 
## sigma = 	1.000	1.071	1.127	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9929	0.0027	0.0044	
## 
## mu = 		 0.000	-3.328	 3.042	
## 
## sigma = 	1.000	1.612	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9964	0.0017	0.0019	
## 
## mu = 		 0.000	-4.213	 3.195	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9456	0.0214	0.0330	
## 
## mu = 		 0.000	-2.996	 2.873	
## 
## sigma = 	1.000	1.000	1.137	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9902	0.0051	0.0047	
## 
## mu = 		 0.000	 2.524	-2.516	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9887	0.0023	0.0090	
## 
## mu = 		 0.000	-4.542	 2.408	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9573	0.0171	0.0256	
## 
## mu = 		 0.000	 2.971	-2.924	
## 
## sigma = 	1.000	1.000	1.266	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9910	0.0071	0.0019	
## 
## mu = 		 0.000	-2.567	 2.257	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9417	0.0307	0.0276	
## 
## mu = 		 0.000	-2.914	 3.100	
## 
## sigma = 	1.00	1.05	1.00	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9936	0.0029	0.0036	
## 
## mu = 		 0.000	 4.636	-4.455	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9936	0.0039	0.0024	
## 
## mu = 		 0.000	 2.695	-3.872	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9187	0.0440	0.0373	
## 
## mu = 		 0.000	 2.790	-2.988	
## 
## sigma = 	1.000	1.171	1.166	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9628	0.0131	0.0241	
## 
## mu = 		 0.000	-3.212	 3.211	
## 
## sigma = 	1.000	1.765	1.667	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9910	0.0034	0.0056	
## 
## mu = 		 0.000	 3.533	-4.147	
## 
## sigma = 	1.000	1.578	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9652	0.0133	0.0215	
## 
## mu = 		 0.000	 3.083	-3.081	
## 
## sigma = 	1.000	1.293	1.527	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9974	0.0016	0.0011	
## 
## mu = 		 0.0000	-3.2379	 0.3495	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9684	0.0165	0.0151	
## 
## mu = 		 0.000	 3.081	-3.444	
## 
## sigma = 	1.000	1.292	1.133	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9973	0.0012	0.0015	
## 
## mu = 		 0.000	-4.514	 1.437	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9348	0.0240	0.0412	
## 
## mu = 		 0.000	-3.007	 2.587	
## 
## sigma = 	1.000	2.266	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9822	0.0069	0.0109	
## 
## mu = 		 0.000	 3.529	-3.241	
## 
## sigma = 	1.000	1.446	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9649	0.0184	0.0167	
## 
## mu = 		 0.000	 3.453	-2.470	
## 
## sigma = 	1.000	1.000	2.387	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.996	0.003	0.001	
## 
## mu = 		 0.000	-2.752	 7.333	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9499	0.0251	0.0249	
## 
## mu = 		 0.000	 3.151	-2.952	
## 
## sigma = 	1.000	1.275	1.172	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9554	0.0238	0.0208	
## 
## mu = 		 0.000	 3.055	-3.009	
## 
## sigma = 	1.000	1.553	1.250	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9529	0.0221	0.0250	
## 
## mu = 		 0.000	 3.005	-2.517	
## 
## sigma = 	1.000	1.278	1.509	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9728	0.0111	0.0161	
## 
## mu = 		 0.000	 3.247	-2.860	
## 
## sigma = 	1.000	1.103	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9717	0.0119	0.0164	
## 
## mu = 		 0.000	-2.676	 3.320	
## 
## sigma = 	1.000	1.000	1.139	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9434	0.0338	0.0228	
## 
## mu = 		 0.000	-2.761	 3.005	
## 
## sigma = 	1.00	1.00	1.27	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9780	0.0128	0.0092	
## 
## mu = 		 0.000	-3.162	 3.198	
## 
## sigma = 	1.000	1.164	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9203	0.0269	0.0528	
## 
## mu = 		 0.000	 3.334	-2.768	
## 
## sigma = 	1.000	1.197	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9466	0.0244	0.0290	
## 
## mu = 		 0.000	 2.922	-3.143	
## 
## sigma = 	1.000	1.093	1.430	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9159	0.0379	0.0463	
## 
## mu = 		 0.000	-2.643	 3.061	
## 
## sigma = 	1.000	1.060	1.331	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9683	0.0125	0.0191	
## 
## mu = 		 0.000	-2.998	 3.052	
## 
## sigma = 	1.000	1.246	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9347	0.0307	0.0346	
## 
## mu = 		 0.000	-3.183	 2.735	
## 
## sigma = 	1.000	1.083	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Using an empirical null with a fitted noiseSD gives a
## substantially different model. Consider rerunning with theonull = FALSE
## and noiseSD = NA.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9081	0.0445	0.0474	
## 
## mu = 		 0.000	 2.848	-2.761	
## 
## sigma = 	1.000	1.024	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9805	0.0085	0.0110	
## 
## mu = 		 0.000	-2.866	 2.683	
## 
## sigma = 	1.000	1.341	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9076	0.0549	0.0375	
## 
## mu = 		 0.000	-2.626	 3.083	
## 
## sigma = 	1.000	1.000	1.341	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9731	0.0132	0.0137	
## 
## mu = 		 0.000	-3.257	 2.686	
## 
## sigma = 	1.000	1.148	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9976	0.0013	0.0011	
## 
## mu = 		 0.000	-2.842	-1.161	
## 
## sigma = 	1.00	1.00	1.32	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9792	0.0120	0.0088	
## 
## mu = 		 0.000	 2.985	-2.922	
## 
## sigma = 	1.000	1.258	1.000	
## 
## noiseSD = 	1	
```

```r
simres1a = basicsim(mixsd, mixpi_alt, niter = 200, nsamp = 10000)
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9456	0.0253	0.0292	
## 
## mu = 		 0.000	 2.987	-3.008	
## 
## sigma = 	1.000	1.214	1.373	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9369	0.0308	0.0323	
## 
## mu = 		 0.000	 2.888	-2.878	
## 
## sigma = 	1.000	1.082	1.190	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9789	0.0105	0.0106	
## 
## mu = 		 0.000	-2.816	 2.476	
## 
## sigma = 	1.000	1.160	1.663	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9659	0.0186	0.0154	
## 
## mu = 		 0.000	 2.774	-2.981	
## 
## sigma = 	1.000	1.304	1.229	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9300	0.0341	0.0360	
## 
## mu = 		 0.000	-3.048	 2.899	
## 
## sigma = 	1.000	1.228	1.047	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9209	0.0405	0.0386	
## 
## mu = 		 0.000	 2.949	-2.835	
## 
## sigma = 	1.000	1.258	1.204	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9861	0.0073	0.0067	
## 
## mu = 		 0.000	 2.128	-2.488	
## 
## sigma = 	1.000	1.996	1.187	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9289	0.0363	0.0348	
## 
## mu = 		 0.000	 2.872	-2.867	
## 
## sigma = 	1.000	1.217	1.072	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9354	0.0309	0.0337	
## 
## mu = 		 0.000	 2.937	-2.991	
## 
## sigma = 	1.000	1.240	1.141	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9854	0.0073	0.0073	
## 
## mu = 		 0.000	-2.791	 2.359	
## 
## sigma = 	1.000	1.129	1.155	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9249	0.0400	0.0352	
## 
## mu = 		 0.000	 2.689	-2.925	
## 
## sigma = 	1.000	1.252	1.164	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9729	0.0135	0.0136	
## 
## mu = 		 0.000	-2.795	 3.115	
## 
## sigma = 	1.000	1.513	1.242	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9661	0.0185	0.0154	
## 
## mu = 		 0.000	 2.886	-2.718	
## 
## sigma = 	1.000	1.415	1.169	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9555	0.0201	0.0244	
## 
## mu = 		 0.000	 2.952	-2.763	
## 
## sigma = 	1.000	1.274	1.201	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9854	0.0074	0.0072	
## 
## mu = 		 0.000	 2.664	-1.569	
## 
## sigma = 	1.000	1.345	2.330	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9873	0.0057	0.0070	
## 
## mu = 		 0.000	-2.536	 2.697	
## 
## sigma = 	1.000	1.563	1.383	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9422	0.0254	0.0324	
## 
## mu = 		 0.000	 3.036	-2.893	
## 
## sigma = 	1.000	1.279	1.128	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9825	0.0089	0.0087	
## 
## mu = 		 0.000	-2.834	 2.986	
## 
## sigma = 	1.000	1.197	1.429	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9398	0.0318	0.0284	
## 
## mu = 		 0.000	-2.913	 3.009	
## 
## sigma = 	1.000	1.107	1.206	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9864	0.0054	0.0081	
## 
## mu = 		 0.000	-2.338	 2.520	
## 
## sigma = 	1.000	1.799	1.337	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9396	0.0287	0.0317	
## 
## mu = 		 0.000	 2.865	-2.770	
## 
## sigma = 	1.000	1.084	1.355	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9868	0.0066	0.0067	
## 
## mu = 		 0.000	-2.393	 2.751	
## 
## sigma = 	1.000	1.510	1.299	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9569	0.0223	0.0209	
## 
## mu = 		 0.000	 3.010	-2.979	
## 
## sigma = 	1.000	1.154	1.123	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9621	0.0203	0.0176	
## 
## mu = 		 0.000	 3.212	-3.092	
## 
## sigma = 	1.000	1.249	1.268	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9510	0.0238	0.0252	
## 
## mu = 		 0.000	 2.927	-2.797	
## 
## sigma = 	1.000	1.248	1.159	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9653	0.0184	0.0163	
## 
## mu = 		 0.000	-2.986	 3.160	
## 
## sigma = 	1.000	1.297	1.244	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9341	0.0321	0.0338	
## 
## mu = 		 0.000	 3.003	-2.614	
## 
## sigma = 	1.000	1.279	1.264	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9793	0.0096	0.0111	
## 
## mu = 		 0.000	-2.858	 2.930	
## 
## sigma = 	1.000	1.198	1.456	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9621	0.0200	0.0179	
## 
## mu = 		 0.000	 2.976	-2.948	
## 
## sigma = 	1.000	1.271	1.250	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9832	0.0083	0.0085	
## 
## mu = 		 0.000	 2.718	-2.696	
## 
## sigma = 	1.000	1.378	1.154	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9610	0.0172	0.0218	
## 
## mu = 		 0.000	-2.942	 2.884	
## 
## sigma = 	1.000	1.208	1.196	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.9885	0.0061	0.0054	
## 
## mu = 		 0.0000	-0.4364	-0.0224	
## 
## sigma = 	1.000	2.194	1.252	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9409	0.0300	0.0291	
## 
## mu = 		 0.000	 2.946	-3.000	
## 
## sigma = 	1.000	1.187	1.067	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.919	0.044	0.037	
## 
## mu = 		 0.000	 2.775	-2.918	
## 
## sigma = 	1.000	1.203	1.300	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9098	0.0473	0.0429	
## 
## mu = 		 0.000	 2.876	-2.945	
## 
## sigma = 	1.000	1.197	1.194	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9465	0.0271	0.0264	
## 
## mu = 		 0.000	 2.844	-2.874	
## 
## sigma = 	1.000	1.290	1.099	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9778	0.0102	0.0120	
## 
## mu = 		 0.000	 3.047	-2.880	
## 
## sigma = 	1.000	1.308	1.409	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9533	0.0236	0.0231	
## 
## mu = 		 0.000	-2.815	 3.017	
## 
## sigma = 	1.000	1.212	1.115	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9898	0.0051	0.0051	
## 
## mu = 		 0.0000	-0.5848	 1.2232	
## 
## sigma = 	1.000	1.023	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9301	0.0352	0.0347	
## 
## mu = 		 0.000	-2.876	 2.874	
## 
## sigma = 	1.000	1.262	1.182	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9373	0.0312	0.0315	
## 
## mu = 		 0.000	 2.922	-2.946	
## 
## sigma = 	1.000	1.105	1.082	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9216	0.0379	0.0405	
## 
## mu = 		 0.000	-2.866	 2.937	
## 
## sigma = 	1.000	1.049	1.178	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9592	0.0184	0.0224	
## 
## mu = 		 0.000	-3.133	 2.870	
## 
## sigma = 	1.000	1.100	1.223	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9375	0.0285	0.0340	
## 
## mu = 		 0.000	 2.931	-2.731	
## 
## sigma = 	1.000	1.108	1.290	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9845	0.0076	0.0078	
## 
## mu = 		 0.000	-2.564	 2.887	
## 
## sigma = 	1.000	1.807	1.454	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9798	0.0096	0.0105	
## 
## mu = 		 0.000	-3.206	 2.676	
## 
## sigma = 	1.000	1.148	1.300	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9887	0.0056	0.0057	
## 
## mu = 		 0.000	-2.266	 2.434	
## 
## sigma = 	1.000	1.235	1.869	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9798	0.0100	0.0101	
## 
## mu = 		 0.000	-2.737	 3.000	
## 
## sigma = 	1.00	1.30	1.31	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9655	0.0199	0.0147	
## 
## mu = 		 0.000	-3.020	 2.891	
## 
## sigma = 	1.000	1.174	1.176	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9586	0.0221	0.0193	
## 
## mu = 		 0.000	-2.752	 2.966	
## 
## sigma = 	1.000	1.441	1.263	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9218	0.0385	0.0397	
## 
## mu = 		 0.000	 2.813	-2.724	
## 
## sigma = 	1.000	1.248	1.312	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9270	0.0363	0.0366	
## 
## mu = 		 0.000	-2.954	 2.895	
## 
## sigma = 	1.000	1.211	1.230	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9687	0.0162	0.0150	
## 
## mu = 		 0.000	-3.037	 2.887	
## 
## sigma = 	1.000	1.131	1.337	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9348	0.0326	0.0326	
## 
## mu = 		 0.000	-2.941	 2.914	
## 
## sigma = 	1.000	1.114	1.207	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	 TRUE	FALSE	
## 
## pi =		0.9899	0.0051	0.0050	
## 
## mu = 		0.0000	0.0999	0.9951	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9770	0.0108	0.0122	
## 
## mu = 		 0.000	-3.091	 2.749	
## 
## sigma = 	1.000	1.221	1.255	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9832	0.0092	0.0076	
## 
## mu = 		 0.000	 2.743	-2.787	
## 
## sigma = 	1.000	1.355	1.236	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9505	0.0255	0.0240	
## 
## mu = 		 0.000	 2.942	-3.062	
## 
## sigma = 	1.000	1.159	1.169	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9486	0.0259	0.0255	
## 
## mu = 		 0.000	 3.002	-2.912	
## 
## sigma = 	1.000	1.237	1.151	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9310	0.0362	0.0328	
## 
## mu = 		 0.000	 2.848	-2.908	
## 
## sigma = 	1.000	1.165	1.157	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9862	0.0076	0.0063	
## 
## mu = 		 0.00	-2.26	 1.13	
## 
## sigma = 	1.000	1.196	2.219	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9761	0.0099	0.0140	
## 
## mu = 		 0.000	 3.246	-2.912	
## 
## sigma = 	1.000	1.142	1.077	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9809	0.0098	0.0093	
## 
## mu = 		 0.000	 2.447	-2.779	
## 
## sigma = 	1.000	1.609	1.232	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9786	0.0103	0.0112	
## 
## mu = 		 0.000	-2.584	 3.007	
## 
## sigma = 	1.000	1.482	1.212	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9235	0.0407	0.0359	
## 
## mu = 		 0.000	-2.692	 2.705	
## 
## sigma = 	1.000	1.273	1.256	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9497	0.0236	0.0267	
## 
## mu = 		 0.000	-2.918	 2.685	
## 
## sigma = 	1.000	1.277	1.259	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9789	0.0101	0.0111	
## 
## mu = 		 0.000	 2.444	-2.780	
## 
## sigma = 	1.000	1.301	1.358	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9405	0.0289	0.0305	
## 
## mu = 		 0.000	 3.078	-2.882	
## 
## sigma = 	1.000	1.119	1.109	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9217	0.0383	0.0400	
## 
## mu = 		 0.000	 2.841	-2.926	
## 
## sigma = 	1.000	1.127	1.115	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9806	0.0094	0.0100	
## 
## mu = 		 0.000	-2.922	 2.695	
## 
## sigma = 	1.000	1.212	1.763	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	 TRUE	FALSE	
## 
## pi =		0.9894	0.0053	0.0053	
## 
## mu = 		 0.0000	-0.1583	-0.3238	
## 
## sigma = 	1.000	1.000	2.074	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9520	0.0272	0.0208	
## 
## mu = 		 0.000	 2.927	-2.965	
## 
## sigma = 	1.000	1.138	1.171	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9880	0.0059	0.0061	
## 
## mu = 		 0.000	-2.536	 2.341	
## 
## sigma = 	1.000	2.000	1.334	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9452	0.0292	0.0255	
## 
## mu = 		 0.000	-2.793	 2.876	
## 
## sigma = 	1.000	1.232	1.157	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9616	0.0215	0.0169	
## 
## mu = 		 0.000	-2.635	 2.937	
## 
## sigma = 	1.000	1.202	1.265	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9205	0.0394	0.0401	
## 
## mu = 		 0.000	 2.769	-2.742	
## 
## sigma = 	1.000	1.146	1.180	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9685	0.0159	0.0156	
## 
## mu = 		 0.000	-2.904	 3.014	
## 
## sigma = 	1.000	1.177	1.474	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9578	0.0210	0.0212	
## 
## mu = 		 0.000	-2.796	 2.716	
## 
## sigma = 	1.000	1.347	1.250	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9550	0.0214	0.0236	
## 
## mu = 		 0.000	-2.873	 2.841	
## 
## sigma = 	1.000	1.159	1.123	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9641	0.0178	0.0181	
## 
## mu = 		 0.000	 2.937	-2.778	
## 
## sigma = 	1.000	1.323	1.389	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9495	0.0244	0.0261	
## 
## mu = 		 0.000	 3.095	-2.949	
## 
## sigma = 	1.000	1.243	1.175	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9266	0.0402	0.0333	
## 
## mu = 		 0.000	-2.906	 2.977	
## 
## sigma = 	1.000	1.269	1.179	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9128	0.0412	0.0460	
## 
## mu = 		 0.000	-2.924	 2.853	
## 
## sigma = 	1.000	1.172	1.159	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9731	0.0129	0.0140	
## 
## mu = 		 0.000	 2.864	-2.807	
## 
## sigma = 	1.000	1.315	1.331	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9375	0.0318	0.0307	
## 
## mu = 		 0.000	-2.915	 2.966	
## 
## sigma = 	1.000	1.245	1.135	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9254	0.0352	0.0394	
## 
## mu = 		 0.000	-2.795	 2.814	
## 
## sigma = 	1.000	1.201	1.232	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9309	0.0346	0.0345	
## 
## mu = 		 0.000	 2.899	-2.868	
## 
## sigma = 	1.000	1.160	1.168	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9815	0.0094	0.0091	
## 
## mu = 		 0.000	 2.463	-3.006	
## 
## sigma = 	1.000	1.496	1.306	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9170	0.0352	0.0479	
## 
## mu = 		 0.000	-2.999	 2.804	
## 
## sigma = 	1.000	1.240	1.204	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9686	0.0155	0.0159	
## 
## mu = 		 0.000	-2.901	 2.990	
## 
## sigma = 	1.000	1.606	1.170	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9381	0.0300	0.0319	
## 
## mu = 		 0.000	-3.034	 2.962	
## 
## sigma = 	1.000	1.145	1.192	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9914	0.0038	0.0048	
## 
## mu = 		 0.0000	 0.5603	-0.8404	
## 
## sigma = 	1.000	2.302	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9740	0.0131	0.0129	
## 
## mu = 		 0.000	-2.892	 2.811	
## 
## sigma = 	1.000	1.438	1.361	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9500	0.0247	0.0253	
## 
## mu = 		 0.000	 2.951	-2.905	
## 
## sigma = 	1.000	1.176	1.309	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9862	0.0063	0.0075	
## 
## mu = 		 0.000	-2.196	 2.511	
## 
## sigma = 	1.000	1.025	1.306	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9530	0.0257	0.0214	
## 
## mu = 		 0.000	-2.710	 2.952	
## 
## sigma = 	1.000	1.297	1.210	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.953	0.023	0.024	
## 
## mu = 		 0.000	-2.951	 3.039	
## 
## sigma = 	1.000	1.241	1.237	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9847	0.0086	0.0067	
## 
## mu = 		 0.000	 2.930	-2.752	
## 
## sigma = 	1.000	1.318	1.215	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9807	0.0104	0.0089	
## 
## mu = 		 0.000	 2.759	-2.802	
## 
## sigma = 	1.000	1.350	1.204	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9677	0.0153	0.0170	
## 
## mu = 		 0.000	 2.920	-2.986	
## 
## sigma = 	1.000	1.209	1.350	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9182	0.0417	0.0401	
## 
## mu = 		 0.000	 2.950	-2.918	
## 
## sigma = 	1.000	1.041	1.124	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9108	0.0414	0.0478	
## 
## mu = 		 0.000	 2.918	-2.906	
## 
## sigma = 	1.000	1.054	1.167	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9447	0.0291	0.0262	
## 
## mu = 		 0.000	-3.066	 2.941	
## 
## sigma = 	1.000	1.161	1.009	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.9883	0.0059	0.0058	
## 
## mu = 		 0.0000	 0.5146	-0.1534	
## 
## sigma = 	1.000	2.091	2.804	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9879	0.0058	0.0064	
## 
## mu = 		 0.000	-1.287	 1.352	
## 
## sigma = 	1.000	1.840	1.727	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9578	0.0207	0.0216	
## 
## mu = 		 0.000	 2.755	-2.985	
## 
## sigma = 	1.000	1.228	1.244	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9330	0.0354	0.0316	
## 
## mu = 		 0.000	 3.019	-2.935	
## 
## sigma = 	1.000	1.167	1.127	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9689	0.0142	0.0169	
## 
## mu = 		 0.000	-2.932	 2.826	
## 
## sigma = 	1.000	1.415	1.196	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9248	0.0386	0.0366	
## 
## mu = 		 0.000	-2.851	 3.008	
## 
## sigma = 	1.000	1.148	1.329	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9864	0.0071	0.0065	
## 
## mu = 		 0.000	-2.458	 2.536	
## 
## sigma = 	1.000	1.525	1.099	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9194	0.0394	0.0412	
## 
## mu = 		 0.000	-2.850	 2.841	
## 
## sigma = 	1.000	1.126	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9756	0.0117	0.0126	
## 
## mu = 		 0.000	-3.083	 2.763	
## 
## sigma = 	1.000	1.189	1.291	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9241	0.0388	0.0371	
## 
## mu = 		 0.000	 2.878	-2.915	
## 
## sigma = 	1.000	1.176	1.260	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9872	0.0061	0.0066	
## 
## mu = 		 0.000	 2.757	-1.620	
## 
## sigma = 	1.000	1.290	1.685	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9233	0.0402	0.0364	
## 
## mu = 		 0.000	 2.945	-2.884	
## 
## sigma = 	1.000	1.130	1.148	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9355	0.0296	0.0348	
## 
## mu = 		 0.000	-2.823	 2.904	
## 
## sigma = 	1.000	1.239	1.232	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9653	0.0164	0.0183	
## 
## mu = 		 0.000	-3.060	 2.973	
## 
## sigma = 	1.000	1.193	1.100	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9720	0.0162	0.0119	
## 
## mu = 		 0.000	 2.820	-2.784	
## 
## sigma = 	1.000	1.526	1.045	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9685	0.0149	0.0166	
## 
## mu = 		 0.000	-2.997	 2.725	
## 
## sigma = 	1.000	1.415	1.296	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9707	0.0138	0.0155	
## 
## mu = 		 0.000	 2.959	-2.802	
## 
## sigma = 	1.000	1.311	1.380	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9316	0.0341	0.0343	
## 
## mu = 		 0.000	 2.912	-2.972	
## 
## sigma = 	1.000	1.107	1.087	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9234	0.0381	0.0385	
## 
## mu = 		 0.000	 2.938	-2.827	
## 
## sigma = 	1.000	1.204	1.014	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9316	0.0336	0.0348	
## 
## mu = 		 0.000	 2.763	-2.829	
## 
## sigma = 	1.000	1.271	1.221	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9479	0.0259	0.0262	
## 
## mu = 		 0.000	 3.074	-3.015	
## 
## sigma = 	1.000	1.147	1.206	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9221	0.0384	0.0395	
## 
## mu = 		 0.000	 2.847	-2.924	
## 
## sigma = 	1.000	1.183	1.185	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9854	0.0084	0.0063	
## 
## mu = 		 0.000	 2.652	-2.248	
## 
## sigma = 	1.000	1.000	1.027	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9428	0.0276	0.0296	
## 
## mu = 		 0.000	 3.029	-3.030	
## 
## sigma = 	1.000	1.144	1.018	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9805	0.0098	0.0097	
## 
## mu = 		 0.000	 2.436	-2.943	
## 
## sigma = 	1.000	1.758	1.110	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9486	0.0257	0.0257	
## 
## mu = 		 0.000	 2.845	-2.933	
## 
## sigma = 	1.000	1.434	1.310	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.978	0.012	0.010	
## 
## mu = 		 0.000	-2.882	 3.073	
## 
## sigma = 	1.000	1.179	1.368	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9777	0.0116	0.0107	
## 
## mu = 		 0.000	-2.976	 3.158	
## 
## sigma = 	1.000	1.331	1.332	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9466	0.0249	0.0285	
## 
## mu = 		 0.000	 3.216	-2.876	
## 
## sigma = 	1.000	1.238	1.091	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9545	0.0226	0.0229	
## 
## mu = 		 0.000	-3.030	 2.989	
## 
## sigma = 	1.000	1.270	1.192	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9741	0.0116	0.0143	
## 
## mu = 		 0.000	-2.694	 2.850	
## 
## sigma = 	1.000	1.502	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9465	0.0266	0.0269	
## 
## mu = 		 0.000	-2.908	 2.672	
## 
## sigma = 	1.000	1.022	1.445	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9650	0.0162	0.0188	
## 
## mu = 		 0.000	-2.831	 2.822	
## 
## sigma = 	1.000	1.334	1.320	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9266	0.0389	0.0345	
## 
## mu = 		 0.000	 2.877	-2.874	
## 
## sigma = 	1.000	1.301	1.198	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9750	0.0112	0.0137	
## 
## mu = 		 0.000	 2.960	-2.739	
## 
## sigma = 	1.000	1.203	1.100	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9587	0.0213	0.0200	
## 
## mu = 		 0.000	 2.855	-3.079	
## 
## sigma = 	1.000	1.268	1.079	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9776	0.0104	0.0120	
## 
## mu = 		 0.000	 2.983	-2.736	
## 
## sigma = 	1.000	1.285	1.231	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9573	0.0214	0.0213	
## 
## mu = 		 0.000	 2.874	-2.649	
## 
## sigma = 	1.000	1.318	1.458	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	 TRUE	FALSE	
## 
## pi =		0.9894	0.0055	0.0051	
## 
## mu = 		0.0000	0.1537	0.4573	
## 
## sigma = 	1.000	2.661	1.309	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9175	0.0396	0.0428	
## 
## mu = 		 0.000	-2.990	 2.828	
## 
## sigma = 	1.000	1.142	1.198	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9865	0.0068	0.0066	
## 
## mu = 		 0.000	-2.544	 2.506	
## 
## sigma = 	1.000	1.797	1.833	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9275	0.0355	0.0370	
## 
## mu = 		 0.000	 3.061	-2.864	
## 
## sigma = 	1.000	1.128	1.355	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9858	0.0069	0.0073	
## 
## mu = 		 0.000	-2.072	 2.341	
## 
## sigma = 	1.000	1.727	1.466	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9544	0.0228	0.0228	
## 
## mu = 		 0.000	 2.775	-2.977	
## 
## sigma = 	1.000	1.228	1.190	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.985	0.007	0.008	
## 
## mu = 		 0.000	-2.252	 2.815	
## 
## sigma = 	1.000	1.202	1.119	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9875	0.0072	0.0053	
## 
## mu = 		 0.000	 2.271	-2.545	
## 
## sigma = 	1.000	1.698	1.470	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9874	0.0060	0.0066	
## 
## mu = 		 0.0000	-0.9812	 1.7210	
## 
## sigma = 	1.000	2.227	1.698	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9800	0.0095	0.0104	
## 
## mu = 		 0.000	-2.762	 2.582	
## 
## sigma = 	1.000	1.439	1.121	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9131	0.0458	0.0411	
## 
## mu = 		 0.000	 2.872	-2.907	
## 
## sigma = 	1.000	1.135	1.150	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9408	0.0305	0.0286	
## 
## mu = 		 0.000	-3.053	 2.968	
## 
## sigma = 	1.000	1.210	1.202	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9821	0.0095	0.0085	
## 
## mu = 		 0.000	-2.447	 2.859	
## 
## sigma = 	1.000	1.123	1.365	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9263	0.0343	0.0394	
## 
## mu = 		 0.000	-2.848	 2.839	
## 
## sigma = 	1.000	1.111	1.141	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9445	0.0270	0.0286	
## 
## mu = 		 0.000	 2.830	-2.947	
## 
## sigma = 	1.000	1.245	1.196	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9598	0.0214	0.0188	
## 
## mu = 		 0.000	-2.689	 2.976	
## 
## sigma = 	1.000	1.304	1.242	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9885	0.0054	0.0061	
## 
## mu = 		 0.0000	-0.4952	 1.6357	
## 
## sigma = 	1.000	2.272	1.474	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9896	0.0053	0.0051	
## 
## mu = 		 0.0000	 0.3547	-1.4221	
## 
## sigma = 	1.000	1.059	1.241	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9852	0.0068	0.0080	
## 
## mu = 		 0.000	 2.926	-2.458	
## 
## sigma = 	1.000	1.382	1.240	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9467	0.0258	0.0275	
## 
## mu = 		 0.00	 2.83	-2.83	
## 
## sigma = 	1.000	1.285	1.257	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9605	0.0213	0.0182	
## 
## mu = 		 0.000	 3.084	-2.829	
## 
## sigma = 	1.000	1.000	1.371	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9356	0.0324	0.0320	
## 
## mu = 		 0.000	 2.736	-2.818	
## 
## sigma = 	1.000	1.239	1.153	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9674	0.0166	0.0160	
## 
## mu = 		 0.000	-2.820	 2.918	
## 
## sigma = 	1.000	1.104	1.338	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9249	0.0369	0.0382	
## 
## mu = 		 0.000	-2.987	 2.954	
## 
## sigma = 	1.000	1.173	1.176	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9689	0.0162	0.0150	
## 
## mu = 		 0.000	 2.716	-2.716	
## 
## sigma = 	1.000	1.423	1.205	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9323	0.0332	0.0346	
## 
## mu = 		 0.000	-2.932	 2.897	
## 
## sigma = 	1.000	1.185	1.125	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9816	0.0096	0.0087	
## 
## mu = 		 0.000	 2.842	-2.393	
## 
## sigma = 	1.000	1.146	1.810	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9876	0.0062	0.0062	
## 
## mu = 		 0.000	-2.165	 2.055	
## 
## sigma = 	1.000	1.357	1.649	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9223	0.0389	0.0387	
## 
## mu = 		 0.000	 2.792	-2.937	
## 
## sigma = 	1.000	1.125	1.185	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9692	0.0151	0.0157	
## 
## mu = 		 0.000	 2.693	-2.995	
## 
## sigma = 	1.000	1.371	1.425	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9332	0.0334	0.0334	
## 
## mu = 		 0.000	-2.929	 2.936	
## 
## sigma = 	1.000	1.206	1.129	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9762	0.0124	0.0113	
## 
## mu = 		 0.000	-2.915	 2.938	
## 
## sigma = 	1.000	1.141	1.326	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9461	0.0253	0.0286	
## 
## mu = 		 0.000	 3.025	-2.958	
## 
## sigma = 	1.000	1.285	1.133	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9157	0.0449	0.0393	
## 
## mu = 		 0.000	-2.897	 2.836	
## 
## sigma = 	1.000	1.190	1.181	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9144	0.0452	0.0404	
## 
## mu = 		 0.000	-2.881	 2.916	
## 
## sigma = 	1.000	1.189	1.128	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9114	0.0452	0.0435	
## 
## mu = 		 0.000	 2.694	-2.737	
## 
## sigma = 	1.000	1.190	1.229	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9712	0.0134	0.0153	
## 
## mu = 		 0.000	 2.929	-2.848	
## 
## sigma = 	1.000	1.201	1.409	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9224	0.0391	0.0386	
## 
## mu = 		 0.000	 2.781	-2.792	
## 
## sigma = 	1.000	1.287	1.230	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9797	0.0098	0.0105	
## 
## mu = 		 0.000	-3.129	 2.853	
## 
## sigma = 	1.000	1.103	1.404	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9855	0.0078	0.0067	
## 
## mu = 		 0.000	-2.361	 3.016	
## 
## sigma = 	1.000	1.259	1.176	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9074	0.0449	0.0476	
## 
## mu = 		 0.000	 2.702	-2.777	
## 
## sigma = 	1.000	1.364	1.213	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9450	0.0297	0.0253	
## 
## mu = 		 0.000	 2.892	-2.965	
## 
## sigma = 	1.000	1.223	1.040	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9716	0.0140	0.0144	
## 
## mu = 		 0.000	 2.923	-3.076	
## 
## sigma = 	1.000	1.228	1.046	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9806	0.0087	0.0107	
## 
## mu = 		 0.000	-2.815	 2.592	
## 
## sigma = 	1.000	1.246	1.708	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9137	0.0419	0.0445	
## 
## mu = 		 0.000	 2.787	-2.957	
## 
## sigma = 	1.000	1.192	1.169	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9565	0.0226	0.0210	
## 
## mu = 		 0.000	 2.968	-3.093	
## 
## sigma = 	1.000	1.210	1.249	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9859	0.0069	0.0072	
## 
## mu = 		 0.000	 2.064	-2.565	
## 
## sigma = 	1.000	1.660	1.024	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9713	0.0145	0.0142	
## 
## mu = 		 0.000	-2.766	 2.807	
## 
## sigma = 	1.000	1.179	1.032	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9748	0.0120	0.0132	
## 
## mu = 		 0.000	-3.271	 2.780	
## 
## sigma = 	1.000	1.000	1.435	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9369	0.0293	0.0338	
## 
## mu = 		 0.000	-2.921	 2.835	
## 
## sigma = 	1.000	1.129	1.306	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9316	0.0356	0.0328	
## 
## mu = 		 0.000	-2.819	 2.985	
## 
## sigma = 	1.000	1.208	1.305	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9217	0.0413	0.0370	
## 
## mu = 		 0.000	-2.791	 2.683	
## 
## sigma = 	1.000	1.113	1.328	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9265	0.0328	0.0407	
## 
## mu = 		 0.000	 2.828	-2.934	
## 
## sigma = 	1.000	1.237	1.079	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9763	0.0112	0.0125	
## 
## mu = 		 0.000	-2.855	 2.890	
## 
## sigma = 	1.00	1.36	1.29	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9122	0.0475	0.0403	
## 
## mu = 		 0.000	 2.687	-2.861	
## 
## sigma = 	1.000	1.189	1.060	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9678	0.0162	0.0160	
## 
## mu = 		 0.000	-3.036	 3.071	
## 
## sigma = 	1.000	1.208	1.235	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9678	0.0175	0.0147	
## 
## mu = 		 0.000	 2.656	-2.973	
## 
## sigma = 	1.000	1.429	1.227	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9411	0.0302	0.0287	
## 
## mu = 		 0.000	 2.919	-2.955	
## 
## sigma = 	1.000	1.151	1.135	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9801	0.0096	0.0103	
## 
## mu = 		 0.000	 3.133	-2.800	
## 
## sigma = 	1.000	1.308	1.309	
## 
## noiseSD = 	1	
```

```r
simres2 = basicsim(c(0, 4), c(1))
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.6207	0.2055	0.1738	
## 
## mu = 		 0.000	-3.959	 4.111	
## 
## sigma = 	1.000	2.590	2.866	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9319	0.0594	0.0087	
## 
## mu = 		 0.000	 1.191	-4.996	
## 
## sigma = 	1.000	5.497	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.6679	0.1574	0.1747	
## 
## mu = 		 0.000	 4.301	-4.192	
## 
## sigma = 	1.000	2.689	2.577	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8613	0.0530	0.0858	
## 
## mu = 		 0.000	 4.752	-4.335	
## 
## sigma = 	1.000	1.978	2.190	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5955	0.2104	0.1941	
## 
## mu = 		 0.000	-3.865	 4.376	
## 
## sigma = 	1.000	2.444	2.598	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7079	0.1637	0.1284	
## 
## mu = 		 0.000	-3.962	 4.031	
## 
## sigma = 	1.000	2.608	2.414	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5424	0.2240	0.2336	
## 
## mu = 		 0.000	-4.058	 4.142	
## 
## sigma = 	1.000	2.741	2.618	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7028	0.1502	0.1470	
## 
## mu = 		 0.000	-4.135	 4.212	
## 
## sigma = 	1.000	2.778	2.672	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9520	0.0207	0.0272	
## 
## mu = 		 0.000	 4.883	-5.104	
## 
## sigma = 	1.000	2.119	1.694	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.6278	0.1993	0.1729	
## 
## mu = 		 0.000	-4.154	 4.380	
## 
## sigma = 	1.000	2.505	2.513	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8082	0.0934	0.0985	
## 
## mu = 		 0.000	-4.538	 4.220	
## 
## sigma = 	1.000	2.400	2.134	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.4921	0.4423	0.0656	
## 
## mu = 		 0.0000	 0.7181	-3.5933	
## 
## sigma = 	1.000	4.822	1.114	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5181	0.2587	0.2233	
## 
## mu = 		 0.000	 4.370	-4.362	
## 
## sigma = 	1.000	2.932	2.449	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9225	0.0400	0.0375	
## 
## mu = 		 0.000	-4.254	 4.468	
## 
## sigma = 	1.000	3.259	2.795	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9619	0.0190	0.0191	
## 
## mu = 		 0.000	 4.931	-4.812	
## 
## sigma = 	1.000	2.354	2.549	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9375	0.0347	0.0277	
## 
## mu = 		 0.000	-4.879	 3.824	
## 
## sigma = 	1.000	2.150	3.142	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9231	0.0366	0.0402	
## 
## mu = 		 0.000	 4.052	-4.710	
## 
## sigma = 	1.000	2.611	1.996	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5765	0.1928	0.2306	
## 
## mu = 		 0.000	 4.223	-3.940	
## 
## sigma = 	1.000	2.467	2.838	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9274	0.0220	0.0506	
## 
## mu = 		 0.000	-4.573	 3.128	
## 
## sigma = 	1.000	2.151	3.716	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5792	0.2222	0.1986	
## 
## mu = 		 0.000	-4.159	 4.422	
## 
## sigma = 	1.000	2.422	2.316	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8495	0.0745	0.0760	
## 
## mu = 		 0.000	-4.152	 4.574	
## 
## sigma = 	1.000	2.552	2.246	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9642	0.0187	0.0171	
## 
## mu = 		 0.000	-2.946	 3.674	
## 
## sigma = 	1.000	1.163	3.483	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.4555	0.2594	0.2851	
## 
## mu = 		 0.000	-4.213	 4.318	
## 
## sigma = 	1.000	2.486	2.430	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5668	0.2060	0.2272	
## 
## mu = 		 0.000	-4.147	 4.239	
## 
## sigma = 	1.000	2.436	2.128	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## null?		 TRUE	 TRUE	FALSE	
## 
## pi =		0.7873	0.1828	0.0299	
## 
## mu = 		 0.0000	-0.0123	 3.2180	
## 
## sigma = 	1.000	4.902	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.6732	0.1597	0.1672	
## 
## mu = 		 0.000	-4.033	 4.266	
## 
## sigma = 	1.000	2.203	2.565	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9173	0.0356	0.0471	
## 
## mu = 		 0.000	-4.742	 4.116	
## 
## sigma = 	1.000	1.976	2.448	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.6459	0.1973	0.1568	
## 
## mu = 		 0.000	 4.319	-4.467	
## 
## sigma = 	1.000	2.425	2.338	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.6592	0.1608	0.1800	
## 
## mu = 		 0.000	-4.595	 4.306	
## 
## sigma = 	1.000	2.345	2.213	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.6992	0.1251	0.1757	
## 
## mu = 		 0.000	 4.609	-3.753	
## 
## sigma = 	1.000	2.391	3.179	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: f(z) misfit = 1.5.  Rerun with increased df
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model
```

```
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.7722	0.0803	0.1476	
## 
## mu = 		 0.000	 3.956	-4.174	
## 
## sigma = 	1.000	2.330	3.073	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7289	0.1214	0.1496	
## 
## mu = 		 0.000	 4.222	-4.525	
## 
## sigma = 	1.000	2.368	2.596	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7244	0.1391	0.1365	
## 
## mu = 		 0.000	 3.957	-4.444	
## 
## sigma = 	1.000	2.334	2.387	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7615	0.1121	0.1263	
## 
## mu = 		 0.000	-4.471	 3.943	
## 
## sigma = 	1.000	2.291	2.256	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5116	0.2439	0.2445	
## 
## mu = 		 0.000	-4.510	 4.342	
## 
## sigma = 	1.000	2.414	2.585	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5049	0.2284	0.2667	
## 
## mu = 		 0.000	 4.369	-4.010	
## 
## sigma = 	1.000	2.466	2.699	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7099	0.1465	0.1437	
## 
## mu = 		 0.000	-3.867	 4.538	
## 
## sigma = 	1.000	2.978	2.491	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5534	0.2338	0.2128	
## 
## mu = 		 0.000	 4.348	-4.186	
## 
## sigma = 	1.000	2.240	2.755	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.5548	0.2338	0.2114	
## 
## mu = 		 0.000	-3.913	 4.298	
## 
## sigma = 	1.000	2.603	2.665	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7694	0.1188	0.1118	
## 
## mu = 		 0.000	 4.250	-4.373	
## 
## sigma = 	1.000	2.252	2.778	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.4463	0.2737	0.2800	
## 
## mu = 		 0.000	 4.080	-4.317	
## 
## sigma = 	1.000	2.372	2.134	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8877	0.0546	0.0577	
## 
## mu = 		 0.000	-4.522	 4.240	
## 
## sigma = 	1.000	2.237	2.728	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9532	0.0167	0.0302	
## 
## mu = 		 0.000	-4.571	 5.059	
## 
## sigma = 	1.000	2.032	1.772	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9719	0.0161	0.0120	
## 
## mu = 		 0.000	-4.200	 4.444	
## 
## sigma = 	1.000	4.657	1.911	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9205	0.0427	0.0368	
## 
## mu = 		 0.000	 3.780	-4.569	
## 
## sigma = 	1.000	1.755	1.798	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9972	0.0014	0.0014	
## 
## mu = 		0.0000	0.4647	0.7202	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.6963	0.1393	0.1644	
## 
## mu = 		 0.000	-4.049	 4.338	
## 
## sigma = 	1.000	2.174	2.200	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.4852	0.2711	0.2438	
## 
## mu = 		 0.000	-4.445	 4.496	
## 
## sigma = 	1.000	2.403	2.413	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.4957	0.2464	0.2578	
## 
## mu = 		 0.000	-4.432	 4.198	
## 
## sigma = 	1.000	2.465	2.414	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9091	0.0208	0.0701	
## 
## mu = 		 0.0000	 3.6510	-0.6115	
## 
## sigma = 	1.000	1.234	4.831	
## 
## noiseSD = 	1	
```

```r
simres3 = basicsim(c(0, 4), c(1), bsd = c(rep(1, 500), rep(10, 500)))
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.7756	0.0667	0.1577	
## 
## mu = 		 0.000	-4.211	 1.754	
## 
## sigma = 	1.000	2.313	4.311	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9589	0.0344	0.0066	
## 
## mu = 		 0.000	 1.359	-4.713	
## 
## sigma = 	1.000	5.271	1.496	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8230	0.0901	0.0869	
## 
## mu = 		 0.000	-4.024	 4.284	
## 
## sigma = 	1.000	2.187	2.904	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9317	0.0406	0.0277	
## 
## mu = 		 0.000	-4.144	 4.423	
## 
## sigma = 	1.000	2.084	2.024	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7654	0.1935	0.0411	
## 
## mu = 		 0.0000	 0.7406	-3.2497	
## 
## sigma = 	1.000	4.937	1.216	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.8366	0.0353	0.1281	
## 
## mu = 		 0.0000	-2.7687	-0.1521	
## 
## sigma = 	1.000	1.000	4.569	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7642	0.1204	0.1155	
## 
## mu = 		 0.000	 4.127	-3.996	
## 
## sigma = 	1.000	2.637	2.600	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8568	0.0751	0.0681	
## 
## mu = 		 0.000	-3.842	 3.989	
## 
## sigma = 	1.000	2.684	2.873	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9801	0.0090	0.0109	
## 
## mu = 		 0.000	 5.483	-4.756	
## 
## sigma = 	1.000	1.475	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8183	0.0935	0.0882	
## 
## mu = 		 0.000	-4.187	 4.272	
## 
## sigma = 	1.000	2.374	2.523	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9126	0.0424	0.0450	
## 
## mu = 		 0.00	-4.79	 4.59	
## 
## sigma = 	1.000	2.129	2.119	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7479	0.1239	0.1282	
## 
## mu = 		 0.000	-3.953	 3.862	
## 
## sigma = 	1.000	2.776	2.582	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7436	0.1283	0.1281	
## 
## mu = 		 0.000	 4.347	-4.353	
## 
## sigma = 	1.000	2.787	2.638	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9568	0.0098	0.0334	
## 
## mu = 		 0.000	 3.213	-1.474	
## 
## sigma = 	1.000	1.000	4.777	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9811	0.0110	0.0080	
## 
## mu = 		 0.000	 4.999	-5.216	
## 
## sigma = 	1.000	2.365	2.597	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9711	0.0218	0.0071	
## 
## mu = 		 0.000	-4.695	 4.427	
## 
## sigma = 	1.000	2.186	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9661	0.0125	0.0213	
## 
## mu = 		 0.000	 4.356	-4.537	
## 
## sigma = 	1.000	1.000	1.873	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.7680	0.0356	0.1964	
## 
## mu = 		 0.0000	-4.0084	 0.3638	
## 
## sigma = 	1.000	1.206	4.900	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	 TRUE	FALSE	
## 
## pi =		0.9668	0.0263	0.0068	
## 
## mu = 		0.0000	0.1122	2.7411	
## 
## sigma = 	1.000	5.152	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8071	0.0977	0.0952	
## 
## mu = 		 0.000	-4.388	 4.771	
## 
## sigma = 	1.000	2.793	2.423	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9220	0.0363	0.0417	
## 
## mu = 		 0.000	-4.205	 4.912	
## 
## sigma = 	1.000	2.422	2.348	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9832	0.0133	0.0035	
## 
## mu = 		 0.000	-2.807	 6.261	
## 
## sigma = 	1.000	1.304	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7154	0.1450	0.1396	
## 
## mu = 		 0.00	-4.04	 4.25	
## 
## sigma = 	1.000	2.615	2.499	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7869	0.0969	0.1162	
## 
## mu = 		 0.000	-4.166	 4.217	
## 
## sigma = 	1.000	2.438	2.271	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9012	0.0787	0.0200	
## 
## mu = 		 0.000	 1.762	-3.527	
## 
## sigma = 	1.000	5.060	1.337	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8259	0.0285	0.1457	
## 
## mu = 		 0.0000	-3.8807	 0.7694	
## 
## sigma = 	1.000	1.000	4.545	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9592	0.0236	0.0171	
## 
## mu = 		 0.000	 4.634	-5.162	
## 
## sigma = 	1.000	2.692	1.481	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: f(z) misfit = 1.7.  Rerun with increased df
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model
```

```
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8225	0.0827	0.0948	
## 
## mu = 		 0.000	-4.410	 4.341	
## 
## sigma = 	1.000	2.608	2.551	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8234	0.0857	0.0909	
## 
## mu = 		 0.000	-4.530	 4.279	
## 
## sigma = 	1.000	2.354	2.271	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8254	0.0323	0.1424	
## 
## mu = 		 0.000	 4.365	-1.379	
## 
## sigma = 	1.000	1.439	4.965	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8935	0.0370	0.0696	
## 
## mu = 		 0.000	 4.104	-4.466	
## 
## sigma = 	1.00	2.61	2.75	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8683	0.0762	0.0555	
## 
## mu = 		 0.000	-4.732	 4.268	
## 
## sigma = 	1.000	2.639	2.334	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8651	0.0647	0.0701	
## 
## mu = 		 0.000	-4.212	 3.893	
## 
## sigma = 	1.000	2.590	2.573	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8579	0.1116	0.0305	
## 
## mu = 		 0.000	-1.305	 4.018	
## 
## sigma = 	1.000	4.367	1.408	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7535	0.1275	0.1190	
## 
## mu = 		 0.000	-4.417	 4.263	
## 
## sigma = 	1.000	2.527	2.549	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7420	0.1267	0.1313	
## 
## mu = 		 0.000	 4.126	-3.902	
## 
## sigma = 	1.000	2.596	2.510	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8658	0.0753	0.0589	
## 
## mu = 		 0.000	-3.369	 4.813	
## 
## sigma = 	1.000	3.173	2.424	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7840	0.1082	0.1079	
## 
## mu = 		 0.000	-4.196	 4.379	
## 
## sigma = 	1.000	2.512	2.001	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7565	0.1242	0.1193	
## 
## mu = 		 0.000	-3.668	 4.151	
## 
## sigma = 	1.000	2.930	2.777	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9032	0.0488	0.0480	
## 
## mu = 		 0.000	 4.346	-4.320	
## 
## sigma = 	1.000	2.223	2.801	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7223	0.1476	0.1300	
## 
## mu = 		 0.000	-4.283	 4.127	
## 
## sigma = 	1.000	2.307	2.488	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9397	0.0294	0.0309	
## 
## mu = 		 0.00	-4.81	 4.62	
## 
## sigma = 	1.000	2.454	2.433	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9760	0.0142	0.0098	
## 
## mu = 		 0.000	 5.253	-4.340	
## 
## sigma = 	1.000	1.745	1.718	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9859	0.0077	0.0064	
## 
## mu = 		 0.000	 4.593	-4.509	
## 
## sigma = 	1.000	2.070	2.731	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.954	0.024	0.022	
## 
## mu = 		 0.000	 3.518	-4.497	
## 
## sigma = 	1.000	1.995	1.964	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9972	0.0014	0.0014	
## 
## mu = 		0.0000	0.4647	0.7202	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## Warning: Null proportion pi0 is small. Consider increasing penalization
## and/or using an empirical null.
```

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.8499	0.0697	0.0804	
## 
## mu = 		 0.000	-4.045	 4.483	
## 
## sigma = 	1.000	2.473	2.161	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: f(z) misfit = 1.6.  Rerun with increased df
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
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
## pi =		0.7497	0.1281	0.1223	
## 
## mu = 		 0.000	-4.411	 4.440	
## 
## sigma = 	1.000	2.256	2.502	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
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
## pi =		0.7265	0.1393	0.1342	
## 
## mu = 		 0.000	-3.966	 4.097	
## 
## sigma = 	1.000	2.605	2.483	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9701	0.0136	0.0163	
## 
## mu = 		 0.000	-2.413	 4.441	
## 
## sigma = 	1.000	5.017	2.006	
## 
## noiseSD = 	1	
```

```r
simres4 = basicsim(c(0, 4), c(1), bsd = c(rep(1, 500), rep(10, 500)), minpi0 = 0.9, 
    seed = 200)
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9949	0.0017	0.0034	
## 
## mu = 		 0.000	 6.049	-3.340	
## 
## sigma = 	1.000	2.870	1.139	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9861	0.0079	0.0060	
## 
## mu = 		 0.000	-3.623	 3.924	
## 
## sigma = 	1.000	1.381	1.854	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9831	0.0102	0.0066	
## 
## mu = 		 0.000	-5.444	 5.692	
## 
## sigma = 	1.000	2.764	2.315	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9958	0.0018	0.0023	
## 
## mu = 		 0.000	-6.904	 2.326	
## 
## sigma = 	1.000	1.218	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9947	0.0031	0.0022	
## 
## mu = 		 0.000	-3.470	 4.144	
## 
## sigma = 	1.000	1.227	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9851	0.0088	0.0061	
## 
## mu = 		 0.000	-6.976	 5.662	
## 
## sigma = 	1.000	1.610	1.072	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9900	0.0036	0.0064	
## 
## mu = 		 0.000	-6.310	 4.703	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9930	0.0029	0.0041	
## 
## mu = 		 0.000	 5.316	-5.070	
## 
## sigma = 	1.000	1.000	1.142	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9840	0.0109	0.0051	
## 
## mu = 		 0.000	 4.392	-6.218	
## 
## sigma = 	1.000	2.772	2.526	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9781	0.0089	0.0130	
## 
## mu = 		 0.000	-4.965	 5.178	
## 
## sigma = 	1.000	2.943	1.922	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9885	0.0048	0.0067	
## 
## mu = 		 0.000	-3.797	 5.908	
## 
## sigma = 	1.000	1.000	2.229	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9961	0.0020	0.0019	
## 
## mu = 		 0.000	-6.604	 7.509	
## 
## sigma = 	1.000	2.481	2.715	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9922	0.0046	0.0032	
## 
## mu = 		 0.000	 4.372	-4.488	
## 
## sigma = 	1.000	1.000	1.008	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9858	0.0092	0.0050	
## 
## mu = 		 0.000	 5.094	-5.443	
## 
## sigma = 	1.000	3.066	1.482	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9902	0.0050	0.0048	
## 
## mu = 		 0.000	 3.697	-3.940	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9796	0.0111	0.0093	
## 
## mu = 		 0.000	 4.144	-5.375	
## 
## sigma = 	1.000	1.772	1.248	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9958	0.0032	0.0010	
## 
## mu = 		 0.000	-6.060	 7.138	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9880	0.0048	0.0072	
## 
## mu = 		 0.000	 5.693	-4.072	
## 
## sigma = 	1.000	2.031	2.255	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9940	0.0038	0.0023	
## 
## mu = 		0.000	4.000	3.818	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: glm.fit: fitted rates numerically 0 occurred
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9808	0.0067	0.0126	
## 
## mu = 		 0.000	-5.838	 2.627	
## 
## sigma = 	1.000	1.690	3.831	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9777	0.0039	0.0183	
## 
## mu = 		 0.000	-3.710	 2.572	
## 
## sigma = 	1.000	1.000	4.878	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9875	0.0048	0.0077	
## 
## mu = 		 0.000	-4.251	 5.391	
## 
## sigma = 	1.000	1.212	2.243	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9854	0.0042	0.0104	
## 
## mu = 		 0.000	 4.955	-3.597	
## 
## sigma = 	1.000	3.219	1.266	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9904	0.0041	0.0055	
## 
## mu = 		 0.000	-5.437	 5.722	
## 
## sigma = 	1.000	2.949	1.427	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9882	0.0034	0.0084	
## 
## mu = 		 0.000	-3.492	 6.949	
## 
## sigma = 	1.000	1.627	2.829	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9725	0.0089	0.0186	
## 
## mu = 		 0.000	 5.931	-4.074	
## 
## sigma = 	1.000	2.356	1.715	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9900	0.0023	0.0077	
## 
## mu = 		 0.000	-4.905	 4.583	
## 
## sigma = 	1.000	1.000	1.246	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9913	0.0060	0.0026	
## 
## mu = 		 0.0000	 0.6215	-3.4060	
## 
## sigma = 	1.000	6.163	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.9970	0.0016	0.0014	
## 
## mu = 		 0.0000	 0.7701	-0.0227	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: longer object length is not a multiple of shorter object length
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9959	0.0023	0.0018	
## 
## mu = 		 0.000	-5.319	 6.127	
## 
## sigma = 	1.000	2.772	1.000	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9820	0.0075	0.0105	
## 
## mu = 		 0.000	 5.854	-3.265	
## 
## sigma = 	1.000	1.454	1.483	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9873	0.0047	0.0079	
## 
## mu = 		 0.000	-7.166	-3.745	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: Assuming known noise noiseSD = 1 . If needed rerun with noiseSD =
## NA to fit noiseSD.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9969	0.0019	0.0012	
## 
## mu = 		 0.000	 4.961	-2.635	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
```

```
## Warning: glm.fit: algorithm did not converge
## Warning: glm.fit: fitted rates numerically 0 occurred
```

```
## Error: system is computationally singular: reciprocal condition number =
## 1.78038e-17
```


Show illustrative example

```r
altcol = "cyan"  #colors to use
nullcol = "blue"
nc = 40  #number of bins in histograms
ncz = 100  # number of bins in z score histograms


plot_FDR_hist = function(sim, iter = 1) {
    hh.pval = sim$pval[[iter]]
    hh.zscore = sim$zscore[[iter]]
    hh.hist = hist(hh.pval, freq = FALSE, xlab = "p value", main = "Distribution of p values", 
        nclass = nc, col = altcol)
    
    hh.q = qvalue(hh.pval)
    abline(h = hh.q$pi0, col = nullcol, lwd = 2)
    
    hh.hist$density = rep(hh.q$pi0, length(hh.hist$density))
    plot(hh.hist, add = TRUE, col = nullcol, freq = FALSE)
    
    abline(v = 0.1, lwd = 2, col = 2)
    
    text(0.05, 1.2, labels = "A", col = 2, cex = 1.2)
    text(0.05, 0.4, labels = "B", col = 2, cex = 1.2)
    text(0.6, 3, labels = paste0("FDR = B/(A+B) =  ", round(hh.q$pi0 * 0.1 * 
        length(hh.pval)/sum(hh.pval < 0.1), 2)), cex = 1.2)
}
plot_FDR_hist(simres1, 1)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 



```r
plot_lfdr_hist = function(sim, iter = 1) {
    require(fdrtool)
    hh.pval = sim$pval[[iter]]
    hh.hist = hist(hh.pval, freq = FALSE, xlab = "p value", main = "Distribution of p values", 
        nclass = nc, col = altcol)
    
    hh.gren = grenander(ecdf(hh.pval))
    abline(h = min(hh.gren$f.knots), col = nullcol, lwd = 2)
    lines(hh.gren$x.knots, hh.gren$f.knots, lwd = 2)
    abline(v = 0.1, lwd = 2, col = 2)
    text(0.1, 0.9, labels = "a", col = 2, cex = 1)
    text(0.1, 0.34, labels = "b", col = 2, cex = 1.2)
    text(0.6, 3, labels = paste0("lfdr = b/(a+b) =  ", round(min(hh.gren$f.knots)/approx(hh.gren$x.knots, 
        hh.gren$f.knots, 0.1)$y, 2)), cex = 1.2)
}
plot_lfdr_hist(simres1, 1)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 




```r
# plot a histogram of z scores, highlighting the alternative distribution of
# z scores that is implied by localfdr values lfdr.
nullalthist = function(z, lfdr, ...) {
    h = hist(z, freq = FALSE, col = nullcol, nclass = ncz, ...)
    avlfdr = unlist(lapply(split(lfdr, cut(z, h$breaks), drop = FALSE), mean))
    h$density = (1 - avlfdr) * h$density
    plot(h, add = TRUE, col = altcol, freq = FALSE)
}

# this one puts the null on the bottom
altnullhist = function(z, lfdr, ...) {
    h = hist(z, freq = FALSE, col = altcol, nclass = ncz, ...)
    avlfdr = unlist(lapply(split(lfdr, cut(z, h$breaks), drop = FALSE), mean))
    h$density = avlfdr * h$density
    plot(h, add = TRUE, col = nullcol, freq = FALSE)
}

plotall_hist = function(sim, iter = 1, histfun = nullalthist) {
    hh.zscore = sim$zscore[[iter]]
    par(mfcol = c(2, 2))
    histfun(hh.zscore, sim$betahat.fdrtool[[iter]]$lfdr, main = "fdrtool")
    histfun(hh.zscore, sim$betahat.locfdr[[iter]]$fdr, main = "locfdr")
    histfun(hh.zscore, sim$betahat.mixfdr[[iter]]$fdr, main = "mixfdr")
    histfun(hh.zscore, sim$betahat.ash.n[[iter]]$lfdr, main = "ash")
    par(mfcol = c(1, 1))
}

# pdf('figures/nullalthist.pdf')
plotall_hist(simres1, 1, nullalthist)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-51.png) 

```r
# dev.off()

# pdf('figures/altnullhist.pdf')
plotall_hist(simres1, 1, altnullhist)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-52.png) 

```r
# dev.off()
```




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

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-61.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-62.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-63.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-64.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-65.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-66.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-67.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-68.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-69.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-610.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-611.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-612.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-613.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-614.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-615.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-616.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-617.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-618.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-619.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-620.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-621.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-622.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-623.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-624.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-625.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-626.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-627.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-628.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-629.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-630.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-631.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-632.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-633.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-634.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-635.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-636.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-637.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-638.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-639.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-640.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-641.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-642.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-643.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-644.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-645.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-646.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-647.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-648.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-649.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-650.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-651.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-652.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-653.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-654.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-655.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-656.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-657.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-658.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-659.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-660.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-661.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-662.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-663.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-664.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-665.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-666.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-667.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-668.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-669.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-670.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-671.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-672.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-673.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-674.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-675.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-676.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-677.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-678.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-679.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-680.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-681.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-682.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-683.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-684.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-685.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-686.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-687.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-688.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-689.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-690.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-691.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-692.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-693.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-694.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-695.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-696.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-697.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-698.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-699.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6100.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6101.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6102.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6103.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6104.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6105.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6106.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6107.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6108.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6109.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6110.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6111.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6112.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6113.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6114.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6115.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6116.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6117.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6118.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6119.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6120.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6121.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6122.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6123.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6124.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6125.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6126.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6127.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6128.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6129.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6130.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6131.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6132.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6133.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6134.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6135.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6136.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6137.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6138.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6139.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6140.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6141.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6142.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6143.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6144.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6145.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6146.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6147.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6148.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6149.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6150.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6151.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6152.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6153.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6154.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6155.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6156.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6157.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6158.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6159.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6160.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6161.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6162.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6163.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6164.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6165.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6166.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6167.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6168.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6169.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6170.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6171.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6172.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6173.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6174.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6175.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6176.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6177.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6178.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6179.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6180.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6181.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6182.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6183.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6184.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6185.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6186.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6187.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6188.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6189.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6190.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6191.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6192.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6193.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6194.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6195.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6196.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6197.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6198.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6199.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6200.png) 

```r
plot_ecdf(simres2)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6201.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6202.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6203.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6204.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6205.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6206.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6207.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6208.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6209.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6210.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6211.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6212.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6213.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6214.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6215.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6216.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6217.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6218.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6219.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6220.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6221.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6222.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6223.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6224.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6225.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6226.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6227.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6228.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6229.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6230.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6231.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6232.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6233.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6234.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6235.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6236.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6237.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6238.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6239.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6240.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6241.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6242.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6243.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6244.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6245.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6246.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6247.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6248.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6249.png) ![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6250.png) 



```r
# Plot pi0 from each method
get_pi0.fdrtool = function(f) {
    f$param[3]
}
get_pi0.locfdr = function(f) {
    f$fp0[1, 3]
}
get_pi0.mixfdr = function(f) {
    f$pi[1]
}

plot_pi0 = function(sims) {
    pi0 = sims$pi0
    pi0_ash.n = unlist(lapply(sims$betahat.ash.n, get_pi0))
    pi0_ash.u = unlist(lapply(sims$betahat.ash.u, get_pi0))
    pi0_fdrtool = unlist(lapply(sims$betahat.fdrtool, get_pi0.fdrtool))
    pi0_locfdr = unlist(lapply(sims$betahat.locfdr, get_pi0.locfdr))
    pi0_mixfdr = unlist(lapply(sims$betahat.mixfdr, get_pi0.mixfdr))
    pi0_qval = unlist(lapply(sims$betahat.qval, "[[", "pi0"))
    
    res = data.frame(pi0 = pi0, qvalue = pi0_qval, mixfdr = pi0_mixfdr, locfdr = pi0_locfdr, 
        fdrtool = pi0_fdrtool, ash.nullbiased = pi0_ash.n, ash.uniform = pi0_ash.u)
    require(reshape2)
    res.melt = melt(res, id.vars = c("pi0"), variable.name = "Method")
    p = ggplot(data = res.melt, aes(pi0, value, colour = Method)) + geom_point(shape = 16) + 
        geom_abline(colour = "black") + xlab("True pi0") + ylab("Estimated pi0")
    print(p + scale_y_continuous(limits = c(0, 1)) + scale_x_continuous(limits = c(0, 
        1)) + coord_equal(ratio = 1))
    
}

plot_pi1 = function(sims) {
    pi0 = sims$pi0
    pi0_ash.n = unlist(lapply(sims$betahat.ash.n, get_pi0))
    pi0_ash.u = unlist(lapply(sims$betahat.ash.u, get_pi0))
    pi0_fdrtool = unlist(lapply(sims$betahat.fdrtool, get_pi0.fdrtool))
    pi0_locfdr = unlist(lapply(sims$betahat.locfdr, get_pi0.locfdr))
    pi0_mixfdr = unlist(lapply(sims$betahat.mixfdr, get_pi0.mixfdr))
    pi0_qval = unlist(lapply(sims$betahat.qval, "[[", "pi0"))
    
    res = data.frame(pi0 = pi0, qvalue = pi0_qval, mixfdr = pi0_mixfdr, locfdr = pi0_locfdr, 
        fdrtool = pi0_fdrtool, ash.nullbiased = pi0_ash.n, ash.uniform = pi0_ash.u)
    require(reshape2)
    res.melt = melt(res, id.vars = c("pi0"), variable.name = "Method")
    p = ggplot(data = res.melt, aes(1 - pi0, log2((1 - value)/(1 - pi0)), colour = Method)) + 
        geom_point(shape = 16) + # geom_abline(colour = 'black') +
    xlab("True pi1") + ylab("log2(Estimated pi1/True pi1)")
    print(p + scale_y_continuous(limits = c(-4, 4)) + scale_x_continuous(limits = c(0, 
        1)))
    
}

pdf("figures/estpi0_sim1.pdf")
plot_pi0(simres1)
```

```
## Loading required package: reshape2
```

```
## Warning: Removed 13 rows containing missing values (geom_point).
```

```r
dev.off()
```

```
## pdf 
##   2
```

```r
pdf("figures/estpi0_sim2.pdf")
plot_pi0(simres2)
dev.off()
```

```
## pdf 
##   2
```


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

```r
plot(betahat.ash.u[[1]]$PosteriorSD, betahat.ash.n[[1]]$PosteriorSD)
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



```r
rmse = function(x, y) {
    sqrt(mean((x - y)^2))
}
get_rmse.ash = function(a, b) {
    rmse(a$PosteriorMean, b)
}
get_rmse.mixfdr = function(a, b) {
    rmse(a$effectSize, b)
}
plot_rmse = function(sims, inczero = FALSE, incbetahat = FALSE) {
    err.ash.n = mapply(get_rmse.ash, sims$betahat.ash.n, sims$beta)
    err.ash.u = mapply(get_rmse.ash, sims$betahat.ash.u, sims$beta)
    err.mixfdr = mapply(get_rmse.mixfdr, sims$betahat.mixfdr, sims$beta)
    err.betahat = mapply(rmse, sims$betahat, sims$beta)
    err.zero = unlist(lapply(sims$beta, rmse, y = 0))
    
    res = data.frame(mixfdr = err.mixfdr, ash.nullbiased = err.ash.n, ash.uniform = err.ash.u)
    if (inczero) {
        res = data.frame(res, zero = err.zero)
    }
    if (incbetahat) {
        res = data.frame(res, betahat = err.betahat)
    }
    require(reshape2)
    res.melt = melt(res, id.vars = c("ash.uniform"), variable.name = "Method")
    
    p = ggplot(data = res.melt, aes(ash.uniform, value, colour = Method)) + 
        geom_point(shape = 16) + geom_abline(colour = "black") + xlab("RMSE (ash.uniform)") + 
        ylab("RMSE (other method)")
    print(p + scale_y_continuous(limits = c(0, max(res))) + scale_x_continuous(limits = c(0, 
        max(res))) + coord_equal(ratio = 1))
}
pdf("figures/rmse_sim1.pdf")
plot_rmse(simres1)
dev.off()
```

```
## pdf 
##   2
```

```r
pdf("figures/rmse_sim2.pdf")
plot_rmse(simres2)
dev.off()
```

```
## pdf 
##   2
```






```r
plot_LR = function(sims) {
    hist(unlist(lapply(sims$betahat.ash.u, get_loglik)) - unlist(lapply(sims$betahat.ash.n, 
        get_loglik)), xlab = "loglik difference", main = "loglik differences for nullbiased prior vs mle", 
        nclass = 10)
}

pdf("figures/logLR.pdf")
plot_LR(simres1)
plot_LR(simres2)
dev.off()
```

```
## pdf 
##   2
```



```r
mean_quant = function(x, mult = 1) {
    x <- na.omit(x)
    sd <- mult * sqrt(var(x))
    mean <- mean(x)
    data.frame(y = median(x), ymin = quantile(x, 0.25), ymax = quantile(x, 0.75))
}

# ptype indicates what type of plot to do maxlfsr controls maximum x axis
# value maxy controls maximum y axis value
plot_lfsr = function(sims, maxlfsr = 0.1, ptype = c("lfsr", "lfsra", "lfdr"), 
    maxy = 1) {
    ptype = match.arg(ptype)
    xlabtype = ifelse(ptype == "lfdr", "lfdr", "lfsr")
    res = list()
    for (i in 1:length(sims)) {
        lfsr.ash.n = unlist(lapply(sims[[i]]$betahat.ash.n, "[[", ptype))
        lfsr.ash.u = unlist(lapply(sims[[i]]$betahat.ash.u, "[[", ptype))
        if (ptype == "lfdr") {
            lfsr.ash.true = unlist(lapply(sims[[i]]$betahat.ash.true, "[[", 
                "lfdr"))
            lfdr.mixfdr = unlist(lapply(sims[[i]]$betahat.mixfdr, "[[", "fdr"))
            lfdr.locfdr = unlist(lapply(sims[[i]]$betahat.locfdr, "[[", "fdr"))
        } else {
            lfsr.ash.true = unlist(lapply(sims[[i]]$betahat.ash.true, "[[", 
                "lfsr"))
        }
        
        subset = lfsr.ash.true < maxlfsr
        
        if (length(subset) > 0) {
            res[[i]] = data.frame(Scenario = i, ash.nullbiased = lfsr.ash.n[subset], 
                ash.uniform = lfsr.ash.u[subset], Bayes = 0.1 * maxlfsr * findInterval(lfsr.ash.true[subset], 
                  seq(0, maxlfsr, length = 11)) - 0.05 * maxlfsr)
            if (ptype == "lfdr") {
                res[[i]] = data.frame(res[[i]], mixfdr = lfdr.mixfdr[subset])
            }
        }
    }
    
    require(reshape2)
    
    res.melt = melt(res, id.vars = c("Bayes", "Scenario"), variable.name = "Method")
    
    
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
        "#D55E00", "#CC79A7")
    labels = c("ash (null-biased)", "ash (uniform)")
    breaks = c("ash.nullbiased", "ash.uniform")
    if (ptype == "lfdr") {
        labels = c("mixfdr", labels)
        breaks = c("mixfdr", breaks)
    }
    
    p = ggplot(data = res.melt, aes(Bayes, value, colour = Method)) + facet_grid(. ~ 
        Scenario) + # scale_fill_manual(values=cbbPalette) +
    # scale_colour_manual(values=cbbPalette) +
    geom_point(size = 1, alpha = 0.1) + stat_smooth(se = FALSE, size = 2) + 
        stat_summary(fun.data = "mean_quant", geom = "ribbon", alpha = 0.25) + 
        geom_abline(colour = "black") + xlab(paste0(xlabtype, " (Bayes)")) + 
        ylab(paste0(ptype, " (Method)"))
    
    print(p + scale_y_continuous(limits = c(0, maxy)) + scale_x_continuous(limits = c(0, 
        maxlfsr)) + scale_colour_manual(values = cbbPalette, breaks = breaks, 
        labels = labels))
}
```



```r
# plots estimated lfdr and lfsr against truth. ptype indicates whether to
# use lfsr or lfsra maxx controls maximum x axis value maxy controls maximum
# y axis value
plot_lfsdr = function(sims, maxx = 0.1, ptype = c("lfsr", "lfsra"), maxy = 1) {
    ptype = match.arg(ptype)
    res = list()
    res.lfsr = list()
    res.lfdr = list()
    
    for (i in 1:length(sims)) {
        lfsr.ash.true = unlist(lapply(sims[[i]]$betahat.ash.true, "[[", "lfsr"))
        lfsr.ash.n = unlist(lapply(sims[[i]]$betahat.ash.n, "[[", ptype))
        lfsr.ash.u = unlist(lapply(sims[[i]]$betahat.ash.u, "[[", ptype))
        
        lfdr.ash.n = unlist(lapply(sims[[i]]$betahat.ash.n, "[[", "lfdr"))
        lfdr.ash.u = unlist(lapply(sims[[i]]$betahat.ash.u, "[[", "lfdr"))
        lfdr.ash.true = unlist(lapply(sims[[i]]$betahat.ash.true, "[[", "lfdr"))
        lfdr.mixfdr = unlist(lapply(sims[[i]]$betahat.mixfdr, "[[", "fdr"))
        lfdr.locfdr = unlist(lapply(sims[[i]]$betahat.locfdr, "[[", "fdr"))
        
        
        subset = lfsr.ash.true < maxx
        
        res.lfsr[[i]] = data.frame(Scenario = i, Measure = "lfsr", ash.nullbiased = lfsr.ash.n[subset], 
            ash.uniform = lfsr.ash.u[subset], Bayes = 0.1 * maxx * findInterval(lfsr.ash.true[subset], 
                seq(0, maxx, length = 11)) - 0.05 * maxx, mixfdr = NA)
        
        subset = lfdr.ash.true < maxx
        res.lfdr[[i]] = data.frame(Scenario = i, Measure = "lfdr", ash.nullbiased = lfdr.ash.n[subset], 
            ash.uniform = lfdr.ash.u[subset], Bayes = 0.1 * maxx * findInterval(lfdr.ash.true[subset], 
                seq(0, maxx, length = 11)) - 0.05 * maxx, mixfdr = lfdr.mixfdr[subset])
        
        
        res[[i]] = rbind(res.lfdr[[i]], res.lfsr[[i]])
    }
    
    require(reshape2)
    
    res.melt = melt(res, id.vars = c("Bayes", "Scenario", "Measure"), variable.name = "Method")
    
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
        "#D55E00", "#CC79A7")
    labels = c("mixfdr", "ash (null-biased)", "ash (uniform)")
    breaks = c("mixfdr", "ash.nullbiased", "ash.uniform")
    
    p = ggplot(data = res.melt, aes(Bayes, value, colour = Method)) + facet_grid(Measure ~ 
        Scenario) + # scale_fill_manual(values=cbbPalette) +
    # scale_colour_manual(values=cbbPalette) +
    geom_point(size = 1, alpha = 0.1) + # stat_smooth(se=FALSE,size=2) +
    stat_summary(fun.data = "mean_quant", geom = "ribbon", alpha = 0.25) + geom_abline(colour = "red", 
        size = 1) + xlab("Truth") + ylab("Estimate")
    
    print(p + scale_y_continuous(limits = c(0, maxy)) + scale_x_continuous(limits = c(0, 
        maxx)) + scale_colour_manual(values = cbbPalette, breaks = breaks, labels = labels))
}
```



```r
png("figures/lfsdr_sim1sim2_blowup.png", height = 427, width = 720)
plot_lfsdr(list(simres1, simres1a, simres2), 0.1, ptype = "lfsra")
```

```
## Warning: Removed 2 rows containing missing values (stat_summary).
## Warning: Removed 4542 rows containing missing values (stat_summary).
## Warning: Removed 50588 rows containing missing values (stat_summary).
## Warning: Removed 14167 rows containing missing values (stat_summary).
## Warning: Removed 2 rows containing missing values (geom_point).
## Warning: Removed 4542 rows containing missing values (geom_point).
## Warning: Removed 50588 rows containing missing values (geom_point).
## Warning: Removed 14167 rows containing missing values (geom_point).
```

```r
dev.off()
```

```
## pdf 
##   2
```




```r
png("figures/lfdr_sim1sim2_blowup.png", height = 160, width = 540)
plot_lfsr(list(simres1, simres1a, simres2), 0.1, ptype = "lfdr")
```

```
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 2 rows containing missing values (stat_smooth).
## Warning: Removed 2 rows containing missing values (stat_summary).
## Warning: Removed 2 rows containing missing values (geom_point).
```

```r
dev.off()
```

```
## pdf 
##   2
```



```r
png("figures/lfsra_sim1sim2_blowup.png", height = 160, width = 540)
plot_lfsr(list(simres1, simres1a, simres2), 0.1, ptype = "lfsra")
```

```
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 1 rows containing missing values (stat_smooth).
```

```
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 1 rows containing missing values (stat_smooth).
```

```
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 292 rows containing missing values (stat_smooth).
## Warning: Removed 307 rows containing missing values (stat_smooth).
## Warning: Removed 1 rows containing missing values (stat_summary).
## Warning: Removed 1 rows containing missing values (stat_summary).
## Warning: Removed 599 rows containing missing values (stat_summary).
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Removed 599 rows containing missing values (geom_point).
```

```r
dev.off()
```

```
## pdf 
##   2
```



```r
png("figures/lfsra_sim1sim2_blowup.png", height = 160, width = 540)
plot_lfsr(list(simres1, simres1a, simres2), 0.1, ptype = "lfsr")
```

```
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is >=1000, so using gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 9 rows containing missing values (stat_smooth).
## Warning: Removed 7 rows containing missing values (stat_smooth).
## Warning: Removed 16 rows containing missing values (stat_summary).
## Warning: Removed 16 rows containing missing values (geom_point).
```

```r
dev.off()
```

```
## pdf 
##   2
```





```r
plot_pi0(simres3)
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-171.png) 

```r
plot_pi1(simres3)
```

```
## Warning: NaNs produced
## Warning: Removed 5 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-172.png) 




QUestion: is the null-biased prior maybe a little too conservative?
Answer: log likelihoods don't suggest they are



```r
# hh.ashtrue = hh.ashz hh.ashtrue$fitted.g$pi =
# c(2/3,1/15,1/15,1/15,1/15,1/15) hh.ashtrue$fitted.g$mean = c(0,0,0,0,0,0)
# hh.ashtrue$fitted.g$sd = sqrt(c(0,1,0.2,0.4,0.8,3))

# loglik(hh.ashtrue,betahat,sebetahat) loglik(hh.ashz,betahat,sebetahat)
```


