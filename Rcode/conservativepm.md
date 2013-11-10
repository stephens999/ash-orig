#look into methods for estimating pi0 conservatively

Do some simulations


```r
basicsim = function(mixsd, mixpi_alt, seedval = 100, nsamp = 1000, niter = 50) {
    require(ashr)
    require(qvalue)
    require(fdrtool)
    require(mixfdr)
    require(locfdr)
    set.seed(seedval)
    beta = list()
    betahatsd = list()
    betahat = list()
    zscore = list()
    pval = list()
    betahat.ash.n = list()
    betahat.ash.u = list()
    betahat.ash.npm = list()
    betahat.qval = list()
    betahat.fdrtool = list()
    betahat.locfdr = list()
    betahat.mixfdr = list()
    pi0 = rep(0, niter)
    for (i in 1:niter) {
        pi0[i] = runif(1, 0, 1)
        mixpi = c(pi0[i], (1 - pi0[i]) * mixpi_alt)
        sd = sample(mixsd, nsamp, prob = mixpi, replace = TRUE)
        beta[[i]] = rnorm(nsamp, 0, sd)
        betahatsd[[i]] = 1
        betahat[[i]] = beta[[i]] + rnorm(nsamp, 0, betahatsd[[i]])
        zscore[[i]] = betahat[[i]]/betahatsd[[i]]
        pval[[i]] = pchisq(zscore[[i]]^2, df = 1, lower.tail = F)
        betahat.ash.n[[i]] = ash(betahat[[i]], betahatsd[[i]], pointmass = TRUE, 
            prior = "nullbiased", gridmult = 2)
        betahat.ash.u[[i]] = ash(betahat[[i]], betahatsd[[i]], pointmass = TRUE, 
            prior = "uniform", gridmult = 2)
        betahat.ash.npm[[i]] = ash(betahat[[i]], betahatsd[[i]], pointmass = FALSE, 
            prior = "uniform", gridmult = 2)
        betahat.qval[[i]] = qvalue(pval[[i]])
        betahat.fdrtool[[i]] = fdrtool(pval[[i]], statistic = "pvalue", plot = FALSE)
        betahat.locfdr[[i]] = locfdr(zscore[[i]], nulltype = 0, plot = 0)
        betahat.mixfdr[[i]] = mixFdr(zscore[[i]], noiseSD = 1, theonull = TRUE, 
            plot = FALSE)
    }
    return(list(beta = beta, betahatsd = betahatsd, betahat = betahat, zscore = zscore, 
        pval = pval, betahat.ash.n = betahat.ash.n, betahat.ash.u = betahat.ash.u, 
        betahat.ash.npm = betahat.ash.npm, betahat.qval = betahat.qval, betahat.fdrtool = betahat.fdrtool, 
        betahat.locfdr = betahat.locfdr, betahat.mixfdr = betahat.mixfdr, pi0 = pi0))
}
mixsd = c(0, 0.25, 0.5, 1, 2)
mixpi_alt = c(0.4, 0.2, 0.2, 0.2)  #mixture proportions under the alternative
simres1 = basicsim(mixsd, mixpi_alt)
```

```
## Loading required package: ashr
## Loading required package: truncnorm
## Loading required package: qvalue
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
## pi =		0.9611	0.0175	0.0215	
## 
## mu = 		 0.000	-3.053	 3.404	
## 
## sigma = 	1.000	1.255	1.720	
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
## pi =		0.9351	0.0417	0.0232	
## 
## mu = 		 0.000	 2.636	-3.357	
## 
## sigma = 	1.000	1.000	1.355	
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
## pi =		0.9545	0.0235	0.0219	
## 
## mu = 		 0.000	-2.683	 3.458	
## 
## sigma = 	1.000	1.000	1.446	
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
## pi =		0.9917	0.0016	0.0067	
## 
## mu = 		 0.000	-1.560	 2.802	
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
## pi =		0.9651	0.0117	0.0232	
## 
## mu = 		 0.000	 3.371	-3.207	
## 
## sigma = 	1.000	1.117	1.455	
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
## pi =		0.9861	0.0086	0.0052	
## 
## mu = 		 0.000	-2.806	 2.687	
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
## pi =		0.944	0.034	0.022	
## 
## mu = 		 0.000	-2.592	 3.549	
## 
## sigma = 	1.000	1.247	1.474	
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
## pi =		0.9262	0.0354	0.0384	
## 
## mu = 		 0.000	 2.868	-2.857	
## 
## sigma = 	1.000	1.478	1.000	
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
## pi =		0.7028	0.1502	0.1470	
## 
## mu = 		 0.000	-4.135	 4.212	
## 
## sigma = 	1.000	2.778	2.672	
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
## pi =		0.9225	0.0400	0.0375	
## 
## mu = 		 0.000	-4.254	 4.468	
## 
## sigma = 	1.000	3.259	2.795	
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
## pi =		0.9619	0.0190	0.0191	
## 
## mu = 		 0.000	 4.931	-4.812	
## 
## sigma = 	1.000	2.354	2.549	
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
## pi =		0.5792	0.2222	0.1986	
## 
## mu = 		 0.000	-4.159	 4.422	
## 
## sigma = 	1.000	2.422	2.316	
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
## pi =		0.7289	0.1214	0.1496	
## 
## mu = 		 0.000	 4.222	-4.525	
## 
## sigma = 	1.000	2.368	2.596	
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
## pi =		0.7694	0.1188	0.1118	
## 
## mu = 		 0.000	 4.250	-4.373	
## 
## sigma = 	1.000	2.252	2.778	
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
## pi =		0.4463	0.2737	0.2800	
## 
## mu = 		 0.000	 4.080	-4.317	
## 
## sigma = 	1.000	2.372	2.134	
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
## pi =		0.9205	0.0427	0.0368	
## 
## mu = 		 0.000	 3.780	-4.569	
## 
## sigma = 	1.000	1.755	1.798	
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


Show illustrative example

```r
plot_FDR_hist = function(sim, iter = 1) {
    hh.pval = sim$pval[[iter]]
    hh.zscore = sim$zscore[[iter]]
    altcol = "cyan"  #colors to use
    nullcol = "blue"
    nc = 40  #number of bins in histograms
    ncz = 100  # number of bins in z score histograms
    
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

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 



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

```
## Error: object 'nc' not found
```




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

```
## Error: object 'ncz' not found
```

```r
# dev.off()

# pdf('figures/altnullhist.pdf')
plotall_hist(simres1, 1, altnullhist)
```

```
## Error: object 'ncz' not found
```

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
    }
}
plot_ecdf(simres1)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-51.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-52.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-53.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-54.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-55.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-56.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-57.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-58.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-59.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-510.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-511.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-512.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-513.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-514.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-515.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-516.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-517.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-518.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-519.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-520.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-521.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-522.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-523.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-524.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-525.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-526.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-527.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-528.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-529.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-530.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-531.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-532.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-533.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-534.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-535.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-536.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-537.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-538.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-539.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-540.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-541.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-542.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-543.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-544.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-545.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-546.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-547.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-548.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-549.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-550.png) 

```r
plot_ecdf(simres2)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-551.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-552.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-553.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-554.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-555.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-556.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-557.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-558.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-559.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-560.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-561.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-562.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-563.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-564.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-565.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-566.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-567.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-568.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-569.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-570.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-571.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-572.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-573.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-574.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-575.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-576.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-577.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-578.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-579.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-580.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-581.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-582.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-583.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-584.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-585.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-586.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-587.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-588.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-589.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-590.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-591.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-592.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-593.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-594.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-595.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-596.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-597.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-598.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-599.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5100.png) 



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

pdf("figures/estpi0_sim1.pdf")
```

```
## Error: cannot open file 'figures/estpi0_sim1.pdf'
```

```r
plot_pi0(simres1)
```

```
## Loading required package: reshape2
```

```
## Error: could not find function "ggplot"
```

```r
dev.off()
```

```
## null device 
##           1
```

```r
pdf("figures/estpi0_sim2.pdf")
```

```
## Error: cannot open file 'figures/estpi0_sim2.pdf'
```

```r
plot_pi0(simres2)
```

```
## Error: could not find function "ggplot"
```

```r
dev.off()
```

```
## Error: cannot shut down device 1 (the null device)
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
```

```
## Error: cannot open file 'figures/rmse_sim1.pdf'
```

```r
plot_rmse(simres1)
```

```
## Error: could not find function "ggplot"
```

```r
dev.off()
```

```
## null device 
##           1
```

```r
pdf("figures/rmse_sim2.pdf")
```

```
## Error: cannot open file 'figures/rmse_sim2.pdf'
```

```r
plot_rmse(simres2)
```

```
## Error: could not find function "ggplot"
```

```r
dev.off()
```

```
## Error: cannot shut down device 1 (the null device)
```






```r
plot_LR = function(sims) {
    hist(unlist(lapply(sims$betahat.ash.u, get_loglik)) - unlist(lapply(sims$betahat.ash.n, 
        get_loglik)), xlab = "loglik difference", main = "loglik differences for nullbiased prior vs mle", 
        nclass = 10)
}

pdf("figures/logLR.pdf")
```

```
## Error: cannot open file 'figures/logLR.pdf'
```

```r
plot_LR(simres1)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-91.png) 

```r
plot_LR(simres2)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-92.png) 

```r
dev.off()
```

```
## null device 
##           1
```



```r
plot(betahat.ash.n[[1]]$lfdr, betahat.ash.u[[1]]$lfdr, xlab = "null-biased prior", 
    ylab = "uniform prior", main = "fdr is not robust to estimated pi0", xlim = c(0, 
        0.2), ylim = c(0, 0.2))
```

```
## Error: object 'betahat.ash.n' not found
```

```r
abline(a = 0, b = 1, col = 2, lwd = 2)
```

```
## Error: plot.new has not been called yet
```



```r
plot(betahat.ash.n[[1]]$lfsr, betahat.ash.u[[1]]$lfsr, xlab = "null-biased prior", 
    ylab = "uniform prior", main = "fsr is more robust to estimated pi0", xlim = c(0, 
        0.2), ylim = c(0, 0.2))
```

```
## Error: object 'betahat.ash.n' not found
```

```r
points(betahat.ash.n[[1]]$lfsr, betahat.ash.npm[[1]]$lfsr, col = 3)
```

```
## Error: object 'betahat.ash.n' not found
```

```r
abline(a = 0, b = 1, col = 2, lwd = 2)
```

```
## Error: plot.new has not been called yet
```



```r
plot(betahat.ash.n[[1]]$lfsra, betahat.ash.u[[1]]$lfsra, xlab = "null-biased prior", 
    ylab = "uniform prior", main = "fsra is still more robust to estimated pi0", 
    xlim = c(0, 0.2), ylim = c(0, 0.2))
```

```
## Error: object 'betahat.ash.n' not found
```

```r
points(betahat.ash.n[[1]]$lfsra, betahat.ash.npm[[1]]$lfsra, col = 3)
```

```
## Error: object 'betahat.ash.n' not found
```

```r
abline(a = 0, b = 1, col = 2, lwd = 2)
```

```
## Error: plot.new has not been called yet
```



```r
plot(betahat.ash.n[[1]]$lfsr, 2 * betahat.ash.npm[[1]]$lfsr, xlab = "null-biased prior", 
    ylab = "uniform prior", main = "fsra is still more robust to estimated pi0", 
    xlim = c(0, 0.2), ylim = c(0, 0.2))
```

```
## Error: object 'betahat.ash.n' not found
```

```r
points(betahat.ash.n[[1]]$lfsra, betahat.ash.npm[[1]]$lfsra, col = 3)
```

```
## Error: object 'betahat.ash.n' not found
```

```r
abline(a = 0, b = 1, col = 2, lwd = 2)
```

```
## Error: plot.new has not been called yet
```


QUestion: is the null-biased prior maybe a little too conservative?
Answer: log likelihoods don't suggest they are


```r
hist(pval[[1]])
```

```
## Error: object 'pval' not found
```





```r
hh.ashtrue = hh.ashz
```

```
## Error: object 'hh.ashz' not found
```

```r
hh.ashtrue$fitted.g$pi = c(2/3, 1/15, 1/15, 1/15, 1/15, 1/15)
```

```
## Error: object 'hh.ashtrue' not found
```

```r
hh.ashtrue$fitted.g$mean = c(0, 0, 0, 0, 0, 0)
```

```
## Error: object 'hh.ashtrue' not found
```

```r
hh.ashtrue$fitted.g$sd = sqrt(c(0, 1, 0.2, 0.4, 0.8, 3))
```

```
## Error: object 'hh.ashtrue' not found
```

```r

loglik(hh.ashtrue, betahat, sebetahat)
```

```
## Error: could not find function "loglik"
```

```r
loglik(hh.ashz, betahat, sebetahat)
```

```
## Error: could not find function "loglik"
```


