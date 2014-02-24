
```r
require("ashr")
```

```
## Loading required package: ashr
## Loading required package: truncnorm
```
  
  ```r
  opts_knit$set(progress = TRUE, verbose = TRUE,root.dir="~/Documents/git/ash/paper/Rcode")
  require(ashr)
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
  #simulate from mixture of normals
  #compare normal, uniform and half-uniform mixtures
  #bsd gives standard deviation of beta
  ashsim=function(mixmean,mixsd,mixpi,bsd=1,seedval = 100,nsamp=1000,niter=50){  
  set.seed(seedval)  
  beta =list()
  betahatsd=list()
  betahat = list()
  fit.ash.n = list()
  fit.ash.u = list()
  fit.ash.hu = list()
  fit.ash.hu.vb=list()
  fit.ash.true = list()
  fit.mixfdr= list()
  fit.ash.fdr.n = list()
  
  
  fit.mixfdr.enull= list()
  fit.mixfdr.J10= list()
  fit.mixfdr.J10P0= list()
  fit.mixfdr.J100= list()
  k = length(mixmean)
  for(i in 1:niter){
    comp = sample(1:k,nsamp,prob=mixpi,replace=TRUE)
    sd = mixsd[comp]
    mean = mixmean[comp]
    beta[[i]] = rnorm(nsamp,mean,sd)
    betahatsd[[i]] = bsd
    betahat[[i]] = beta[[i]]+rnorm(nsamp,0,betahatsd[[i]])
    fit.ash.n[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="normal",method="shrink")
    fit.ash.u[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="uniform", method="shrink")
    fit.ash.hu[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="halfuniform", method="shrink")
    fit.ash.hu.vb[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="halfuniform",method="shrink",VB=TRUE)
    fit.ash.true[[i]] = ash(betahat[[i]],betahatsd[[i]],g=normalmix(mixpi,mixmean,mixsd))
    fit.mixfdr[[i]] = mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=TRUE,plot=FALSE)
    
    fit.ash.fdr.n[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="normal",method="fdr")
    #fit.ash.fdr.u[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="uniform", method="fdr")
    #fit.ash.fdr.hu[[i]] = ash(betahat[[i]],betahatsd[[i]],mixcompdist="halfuniform", method="fdr")
    fit.mixfdr.enull[[i]] = mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=FALSE,plot=FALSE)
    fit.mixfdr.J10[[i]] = mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=TRUE,plot=FALSE,J=10)
     #fit.mixfdr.J100[[i]] = mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=TRUE,plot=FALSE,J=100)
    fit.mixfdr.J10P0[[i]] = try(mixFdr(betahat[[i]]/betahatsd[[i]],noiseSD=1,theonull=TRUE,plot=FALSE,J=10,P=0))
  }
  return(list(beta =beta,
  betahatsd=betahatsd,
  betahat = betahat,
  fit.ash.n = fit.ash.n,
  fit.ash.u = fit.ash.u,
  fit.ash.hu = fit.ash.hu,
  fit.ash.hu.vb=fit.ash.hu.vb,
  fit.mixfdr = fit.mixfdr,
  fit.ash.fdr.n = fit.ash.fdr.n,
  fit.ash.true=fit.ash.true,
   # fit.ash.fdr.u = fit.ash.fdr.u,
   # fit.ash.fdr.hu = fit.ash.fdr.hu,
  fit.mixfdr.enull = fit.mixfdr.enull,
  fit.mixfdr.J10 = fit.mixfdr.J10,
  #fit.mixfdr.J100= fit.mixfdr.J100,
  fit.mixfdr.J10P0=fit.mixfdr.J10P0))
  }    
  ```



```r
sim1 = ashsim(c(0, 0, 0), c(1, 1, 2), c(1/3, 1/3, 1/3), niter = 20, nsamp = 1000)
```

```
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8024	0.0925	0.1051	
## 
## mu = 		 0.00	-2.67	 2.76	
## 
## sigma = 	1.000	1.206	1.152	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9829	0.0099	0.0072	
## 
## mu = 		 0.0465	 4.6458	-5.0453	
## 
## sigma = 	1.622	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7871	0.0555	0.0215	0.0091	0.0043	0.0219	0.0273	0.0228	0.0240	0.0265	
## 
## mu = 		 0.000	 2.409	-2.394	-4.932	 3.563	-2.394	 2.407	-2.394	-2.394	 3.375	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.502	1.000	1.000	1.000	1.000	1.347	
## 
## noiseSD = 	1	
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
## pi =		0.5636	0.0427	0.0253	0.0254	0.0367	0.1615	0.0138	0.0497	0.0360	0.0452	
## 
## mu = 		 0.000	 3.485	 1.714	 1.715	 1.750	-1.690	 1.681	-2.820	 1.747	 1.801	
## 
## sigma = 	1.000	1.321	1.000	1.000	1.000	1.000	1.000	1.687	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8258	0.0849	0.0892	
## 
## mu = 		 0.000	 2.796	-2.622	
## 
## sigma = 	1.000	1.247	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.9874	0.0111	0.0015	
## 
## mu = 		-0.0531	 4.7298	-0.1548	
## 
## sigma = 	1.600	1.000	1.554	
## 
## noiseSD = 	1	
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
## pi =		0.8055	0.0082	0.0282	0.0635	0.0039	0.0377	0.0121	0.0112	0.0093	0.0203	
## 
## mu = 		 0.000	 3.560	-2.857	 2.229	 3.014	-2.376	 3.418	-2.475	 3.695	-2.373	
## 
## sigma = 	1.000	1.360	1.099	1.000	1.355	1.000	1.352	1.000	1.352	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5273	0.0182	0.0613	0.0318	0.0037	0.0500	0.0270	0.0763	0.1375	0.0669	
## 
## mu = 		 0.000	-1.460	-2.607	-1.474	 1.227	-1.488	 4.282	-1.511	 1.541	 1.538	
## 
## sigma = 	1.000	1.000	1.270	1.000	1.206	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7901	0.1108	0.0990	
## 
## mu = 		 0.000	-2.687	 2.702	
## 
## sigma = 	1.000	1.226	1.126	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9897	0.0060	0.0043	
## 
## mu = 		-0.0543	-5.4553	 5.2143	
## 
## sigma = 	1.683	1.000	1.086	
## 
## noiseSD = 	1	
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
## pi =		0.7721	0.0402	0.0305	0.0406	0.0425	0.0075	0.0508	0.0035	0.0095	0.0028	
## 
## mu = 		 0.000	 2.490	-2.392	-2.392	-2.392	 4.157	 2.490	 3.400	-5.130	 3.353	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.526	1.000	1.366	1.000	1.468	
## 
## noiseSD = 	1	
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
## pi =		0.5604	0.0638	0.0735	0.0082	0.0554	0.0030	0.0116	0.1275	0.0565	0.0399	
## 
## mu = 		 0.0000	 2.4318	-1.7862	 1.8526	-2.8610	-0.8845	-1.7060	 1.9317	-1.7637	-1.7487	
## 
## sigma = 	1.000	1.681	1.000	1.000	1.757	1.957	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8041	0.0932	0.1027	
## 
## mu = 		 0.000	-2.714	 2.766	
## 
## sigma = 	1.000	1.202	1.073	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9931	0.0057	0.0012	
## 
## mu = 		 0.0777	-5.2231	-0.5525	
## 
## sigma = 	1.687	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7937	0.0049	0.0043	0.0536	0.0055	0.0140	0.0172	0.0604	0.0084	0.0381	
## 
## mu = 		 0.000	-4.269	 2.611	 2.500	 3.257	-2.420	-3.442	-2.420	 4.492	 2.493	
## 
## sigma = 	1.000	1.426	1.000	1.000	1.165	1.000	1.431	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## null?		 TRUE	FALSE	FALSE	 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	
## 
## pi =		0.6109	0.0287	0.0092	0.0032	0.0442	0.0801	0.0038	0.0161	0.0205	0.1834	
## 
## mu = 		 0.0000	-1.8620	-1.8299	-0.1585	-2.9314	-1.9232	-1.6619	-1.8499	 1.7149	 2.1330	
## 
## sigma = 	1.000	1.000	1.000	2.372	1.703	1.000	1.000	1.000	1.000	1.324	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7878	0.1130	0.0992	
## 
## mu = 		 0.000	-2.603	 2.706	
## 
## sigma = 	1.000	1.179	1.364	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9893	0.0059	0.0048	
## 
## mu = 		-0.0427	-5.2404	 5.9776	
## 
## sigma = 	1.673	1.000	2.106	
## 
## noiseSD = 	1	
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
## pi =		0.7683	0.0559	0.0046	0.0068	0.0199	0.0309	0.0039	0.0590	0.0040	0.0466	
## 
## mu = 		 0.000	-2.288	 4.352	-4.959	 2.509	 2.509	 6.048	-2.288	-4.786	 2.509	
## 
## sigma = 	1.000	1.000	1.277	1.000	1.000	1.000	2.342	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5549	0.0364	0.0373	0.0285	0.0151	0.0420	0.2247	0.0256	0.0178	0.0177	
## 
## mu = 		 0.000	 2.018	 2.022	 1.988	 1.942	 2.076	-1.760	 2.857	-4.628	 1.953	
## 
## sigma = 	1.000	1.012	1.014	1.004	1.000	1.033	1.000	2.756	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8241	0.0828	0.0931	
## 
## mu = 		 0.000	-2.907	 2.830	
## 
## sigma = 	1.000	1.030	1.101	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.9978	0.0010	0.0013	
## 
## mu = 		 0.0268	 7.5066	-0.0062	
## 
## sigma = 	1.702	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8133	0.0031	0.0825	0.0022	0.0167	0.0335	0.0054	0.0018	0.0199	0.0216	
## 
## mu = 		 0.000	 3.023	-2.720	 4.957	 2.720	 2.744	-4.512	 1.500	 2.724	 2.726	
## 
## sigma = 	1.000	1.050	1.000	2.169	1.000	1.000	1.297	3.041	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6742	0.0030	0.0431	0.0065	0.0193	0.0225	0.1234	0.0803	0.0029	0.0248	
## 
## mu = 		 0.0000	-1.3989	-2.2126	 2.8726	 1.8895	-3.2167	 2.3702	-2.2798	-0.7915	 1.8995	
## 
## sigma = 	1.000	1.905	1.014	2.822	1.000	1.543	1.186	1.019	2.094	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8459	0.0769	0.0773	
## 
## mu = 		 0.000	-2.749	 2.801	
## 
## sigma = 	1.000	1.055	1.277	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9907	0.0012	0.0081	
## 
## mu = 		-0.0192	 0.2822	 4.9384	
## 
## sigma = 	1.579	1.000	1.107	
## 
## noiseSD = 	1	
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
## pi =		0.8328	0.0470	0.0074	0.0134	0.0215	0.0035	0.0309	0.0099	0.0039	0.0297	
## 
## mu = 		 0.000	 2.430	-3.924	-2.575	 2.425	 3.434	-2.577	 4.286	 3.449	-2.576	
## 
## sigma = 	1.000	1.000	1.218	1.000	1.000	1.451	1.000	1.418	1.371	1.000	
## 
## noiseSD = 	1	
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
##  [1] 0.9120258 0.0003671 0.0017498 0.0000388 0.0163380 0.0000000 0.0300575
##  [8] 0.0011149 0.0002140 0.0380942
##  [1]  0.06627  2.67122  4.84078  0.37851 -2.84065  0.04496  3.59582
##  [8]  4.16505  2.14644 -2.84129
##  [1] 1.352 2.215 1.648 1.279 1.279 1.379 1.279 1.830 2.334 1.279
##           [,1]      [,2]  [,3]
##  [1,] 0.990000  0.000000 1.279
##  [2,] 0.001111  0.836638 1.279
##  [3,] 0.001111  1.199349 1.279
##  [4,] 0.001111  0.285839 1.279
##  [5,] 0.001111 -1.606012 1.279
##  [6,] 0.001111 -0.003665 1.279
##  [7,] 0.001111  3.371097 1.279
##  [8,] 0.001111  1.104524 1.279
##  [9,] 0.001111  0.734988 1.279
## [10,] 0.001111 -2.122183 1.279
##       [,1] [,2]
##  [1,] -Inf  Inf
##  [2,] -Inf  Inf
##  [3,] -Inf  Inf
##  [4,] -Inf  Inf
##  [5,] -Inf  Inf
##  [6,] -Inf  Inf
##  [7,] -Inf  Inf
##  [8,] -Inf  Inf
##  [9,] -Inf  Inf
## [10,] -Inf  Inf
##        [,1] [,2]
##  [1,] 1.279  Inf
##  [2,] 1.279  Inf
##  [3,] 1.279  Inf
##  [4,] 1.279  Inf
##  [5,] 1.279  Inf
##  [6,] 1.279  Inf
##  [7,] 1.279  Inf
##  [8,] 1.279  Inf
##  [9,] 1.279  Inf
## [10,] 1.279  Inf
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7975	0.1007	0.1017	
## 
## mu = 		 0.000	-2.648	 2.687	
## 
## sigma = 	1.000	1.071	1.062	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9951	0.0023	0.0026	
## 
## mu = 		 0.0167	 6.1074	-5.7252	
## 
## sigma = 	1.659	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7888	0.0032	0.0318	0.0027	0.0168	0.0388	0.0027	0.0517	0.0329	0.0307	
## 
## mu = 		 0.000	-5.644	-2.528	 2.718	 2.556	-2.528	 6.057	 2.556	 2.556	-2.528	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.6370	0.0684	0.0075	0.0225	0.0143	0.0102	0.0236	0.0404	0.0830	0.0932	
## 
## mu = 		 0.000	 2.069	 4.532	-1.992	-1.973	-4.318	 2.011	-2.010	 2.246	-2.177	
## 
## sigma = 	1.000	1.015	1.620	1.000	1.000	1.410	1.015	1.000	1.009	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7935	0.0998	0.1067	
## 
## mu = 		 0.000	-2.658	 2.766	
## 
## sigma = 	1.000	1.179	1.214	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9884	0.0058	0.0058	
## 
## mu = 		 0.0491	-5.4306	 5.5351	
## 
## sigma = 	1.661	1.000	1.221	
## 
## noiseSD = 	1	
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
## pi =		0.7768	0.0510	0.0493	0.0031	0.0031	0.0027	0.0373	0.0071	0.0066	0.0628	
## 
## mu = 		 0.000	-2.411	-2.411	 3.610	 3.606	 3.510	 2.498	-5.352	 4.710	 2.498	
## 
## sigma = 	1.000	1.000	1.000	1.576	1.576	1.639	1.000	1.000	1.654	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5975	0.0165	0.1327	0.0339	0.1185	0.0160	0.0333	0.0087	0.0040	0.0391	
## 
## mu = 		 0.0000	 1.7798	 2.1188	-1.9263	-1.9345	 4.4323	-3.0669	-1.8887	 0.7215	 1.8192	
## 
## sigma = 	1.000	1.000	1.043	1.000	1.000	1.533	1.841	1.000	2.565	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7988	0.1019	0.0993	
## 
## mu = 		 0.000	-2.621	 2.734	
## 
## sigma = 	1.000	1.155	1.309	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9874	0.0063	0.0062	
## 
## mu = 		 0.0277	-5.1778	 5.7528	
## 
## sigma = 	1.639	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7799	0.0054	0.0228	0.0329	0.0070	0.0258	0.0749	0.0082	0.0048	0.0382	
## 
## mu = 		 0.000	-2.383	-2.346	 2.462	 5.661	 2.462	-2.346	-5.071	 2.913	 2.462	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.019	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5436	0.0825	0.0032	0.0416	0.0809	0.1276	0.0228	0.0214	0.0546	0.0218	
## 
## mu = 		 0.000	-1.772	 0.280	-2.865	 1.708	 1.752	 4.371	-1.759	-1.764	-1.759	
## 
## sigma = 	1.000	1.000	2.204	1.783	1.000	1.011	1.350	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7828	0.1029	0.1143	
## 
## mu = 		 0.000	-2.647	 2.525	
## 
## sigma = 	1.00	1.09	1.00	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9947	0.0013	0.0040	
## 
## mu = 		 0.0046	 3.2520	-4.3966	
## 
## sigma = 	1.700	1.855	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7654	0.0375	0.0122	0.0438	0.0037	0.0068	0.0517	0.0161	0.0412	0.0215	
## 
## mu = 		 0.000	-2.222	 2.410	 2.416	 3.834	-2.217	-2.924	-2.204	 2.414	 2.412	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.641	1.000	1.180	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4600	0.0257	0.0624	0.0822	0.0123	0.1260	0.0821	0.0699	0.0751	0.0044	
## 
## mu = 		 0.000	 1.835	-1.229	 1.876	 3.353	 1.904	-2.551	-1.230	-1.471	-1.570	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.453	1.000	1.418	1.000	1.000	1.518	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8309	0.0816	0.0876	
## 
## mu = 		 0.000	 2.739	-2.704	
## 
## sigma = 	1.000	1.000	1.149	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9943	0.0019	0.0039	
## 
## mu = 		 0.0077	-0.3259	-5.7438	
## 
## sigma = 	1.612	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8187	0.0040	0.0129	0.0568	0.0038	0.0213	0.0092	0.0288	0.0291	0.0153	
## 
## mu = 		 0.000	-5.770	 2.684	-2.501	-2.719	-2.501	-2.501	 2.688	 2.688	 2.686	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.6459	0.0342	0.0050	0.0033	0.0362	0.1327	0.0450	0.0449	0.0145	0.0384	
## 
## mu = 		 0.000	-2.490	-1.796	-2.123	 1.753	 2.208	-1.992	-1.921	-1.897	-1.918	
## 
## sigma = 	1.000	2.044	1.000	2.247	1.000	1.181	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7822	0.1187	0.0991	
## 
## mu = 		 0.000	-2.720	 2.825	
## 
## sigma = 	1.000	1.111	1.197	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9971	0.0010	0.0018	
## 
## mu = 		-0.0344	-3.3161	-6.0278	
## 
## sigma = 	1.796	1.985	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7632	0.0337	0.0652	0.0180	0.0090	0.0020	0.0034	0.0230	0.0177	0.0648	
## 
## mu = 		 0.0000	-2.5150	-2.5156	-2.5147	-4.4949	-0.8096	 3.3018	 3.7735	 2.4097	 2.4097	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.392	2.964	1.382	1.361	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5271	0.0585	0.0064	0.0065	0.0584	0.1770	0.0229	0.0526	0.0704	0.0203	
## 
## mu = 		 0.000	-1.825	-1.851	-1.247	-2.624	 1.717	 3.662	-1.818	-1.842	 3.448	
## 
## sigma = 	1.000	1.000	2.267	2.433	1.627	1.000	1.320	1.000	1.001	1.422	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8168	0.0887	0.0945	
## 
## mu = 		 0.000	-2.628	 2.770	
## 
## sigma = 	1.000	1.075	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9928	0.0036	0.0036	
## 
## mu = 		 0.0196	-5.0532	 4.3851	
## 
## sigma = 	1.631	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8003	0.0338	0.0168	0.0052	0.0637	0.0171	0.0331	0.0036	0.0129	0.0134	
## 
## mu = 		 0.000	-2.299	 2.636	 3.689	 2.676	-2.299	-2.299	-3.141	-3.584	 2.632	
## 
## sigma = 	1.000	1.000	1.000	1.183	1.000	1.000	1.000	1.347	1.398	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5887	0.0058	0.0515	0.0506	0.0686	0.0301	0.1260	0.0074	0.0468	0.0245	
## 
## mu = 		 0.000	 1.680	 2.772	 2.005	-2.166	 1.845	-1.667	 1.710	 1.933	-1.635	
## 
## sigma = 	1.000	1.000	1.321	1.005	1.659	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8067	0.0937	0.0997	
## 
## mu = 		 0.000	-2.654	 2.565	
## 
## sigma = 	1.000	1.000	1.076	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9971	0.0013	0.0015	
## 
## mu = 		 0.0277	-3.7263	 6.1938	
## 
## sigma = 	1.643	1.184	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7915	0.0027	0.0059	0.0461	0.0163	0.0173	0.0360	0.0159	0.0592	0.0090	
## 
## mu = 		 0.000	-3.406	 4.435	-2.583	-2.583	-2.583	 2.363	-2.583	 2.366	 2.364	
## 
## sigma = 	1.000	1.538	1.522	1.000	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6176	0.0426	0.0202	0.0255	0.0188	0.0208	0.1218	0.0683	0.0058	0.0585	
## 
## mu = 		 0.000	-2.172	-2.113	-2.132	 3.388	-2.116	 1.881	 1.827	-2.927	-2.207	
## 
## sigma = 	1.000	1.028	1.024	1.026	1.659	1.025	1.000	1.000	1.792	1.029	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8056	0.0966	0.0979	
## 
## mu = 		 0.000	-2.872	 2.960	
## 
## sigma = 	1.000	1.161	1.138	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9867	0.0078	0.0056	
## 
## mu = 		-0.0126	 4.8967	-5.0088	
## 
## sigma = 	1.699	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7997	0.0448	0.0666	0.0196	0.0061	0.0144	0.0202	0.0041	0.0131	0.0115	
## 
## mu = 		 0.000	 2.722	-2.616	 2.718	 2.954	 4.056	-2.558	 3.268	-4.447	 2.720	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.327	1.000	1.211	1.045	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6457	0.0862	0.0025	0.0045	0.0856	0.0272	0.0489	0.0428	0.0023	0.0542	
## 
## mu = 		 0.0000	 2.1220	-1.3475	-1.6591	 2.7439	-1.8746	-1.9847	-1.9454	-0.5457	-3.1761	
## 
## sigma = 	1.000	1.000	1.631	1.000	1.485	1.000	1.000	1.000	1.996	1.411	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8042	0.1088	0.0871	
## 
## mu = 		 0.000	 2.715	-2.824	
## 
## sigma = 	1.000	1.062	1.085	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.9966	0.0022	0.0012	
## 
## mu = 		0.0379	5.6665	0.0127	
## 
## sigma = 	1.726	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7862	0.0091	0.0198	0.0198	0.0037	0.0116	0.1042	0.0110	0.0188	0.0159	
## 
## mu = 		 0.000	 2.576	-2.502	-3.589	 5.454	-2.475	 2.550	-2.474	-2.483	-2.478	
## 
## sigma = 	1.000	1.000	1.000	1.177	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5901	0.0821	0.0459	0.0497	0.1103	0.0512	0.0079	0.0124	0.0315	0.0188	
## 
## mu = 		 0.000	-1.626	 1.995	 2.340	-2.275	 2.053	-2.327	 3.998	 1.926	 1.879	
## 
## sigma = 	1.000	1.000	1.000	1.125	1.395	1.002	1.745	1.468	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8040	0.0931	0.1029	
## 
## mu = 		 0.000	-2.734	 2.596	
## 
## sigma = 	1.000	1.163	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9933	0.0011	0.0056	
## 
## mu = 		 0.0565	-0.2088	-5.3175	
## 
## sigma = 	1.641	1.591	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8014	0.0217	0.0310	0.0114	0.0036	0.0135	0.0898	0.0053	0.0167	0.0056	
## 
## mu = 		 0.000	-2.566	-2.566	-2.566	-2.772	 2.593	 2.593	-2.576	-2.566	-5.428	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.6332	0.0176	0.0793	0.0107	0.0643	0.0713	0.0962	0.0096	0.0145	0.0032	
## 
## mu = 		 0.0000	-1.9569	 2.1527	-1.9342	-2.4360	-1.9913	 2.1739	 1.7618	-1.9496	 0.8758	
## 
## sigma = 	1.000	1.000	1.088	1.000	1.719	1.000	1.095	1.000	1.000	1.575	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8165	0.1027	0.0808	
## 
## mu = 		 0.000	 2.643	-2.841	
## 
## sigma = 	1.000	1.009	1.161	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9889	0.0075	0.0036	
## 
## mu = 		 0.0884	-4.9627	 4.1721	
## 
## sigma = 	1.623	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7962	0.0329	0.0179	0.0043	0.0375	0.0231	0.0209	0.0315	0.0022	0.0335	
## 
## mu = 		 0.0000	-2.4366	 2.3318	 2.7352	 2.8750	-3.5190	 2.3337	 2.3422	-0.3258	-2.4366	
## 
## sigma = 	1.000	1.000	1.000	1.126	1.115	1.415	1.000	1.000	2.645	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5648	0.1117	0.0081	0.1227	0.0112	0.0061	0.0312	0.0128	0.0524	0.0790	
## 
## mu = 		 0.000	-1.730	-2.160	 2.187	 1.534	 1.592	 1.563	 1.539	 1.578	-2.245	
## 
## sigma = 	1.000	1.000	1.894	1.332	1.000	1.000	1.000	1.000	1.000	1.661	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7868	0.1145	0.0987	
## 
## mu = 		 0.000	-2.683	 2.589	
## 
## sigma = 	1.000	1.364	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9895	0.0027	0.0078	
## 
## mu = 		-0.0633	 4.3841	-5.6078	
## 
## sigma = 	1.668	1.000	1.329	
## 
## noiseSD = 	1	
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
## pi =		0.7715	0.0175	0.0390	0.0085	0.0062	0.0144	0.0659	0.0609	0.0060	0.0100	
## 
## mu = 		 0.000	 2.445	-2.336	 2.453	 3.973	-3.764	-2.336	 2.451	-5.337	 2.444	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.073	1.454	1.000	1.000	1.672	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5148	0.0813	0.0806	0.0032	0.0064	0.0459	0.0995	0.0711	0.0062	0.0909	
## 
## mu = 		 0.000	-2.094	-1.615	 1.058	-3.723	 1.762	-1.617	 2.295	-3.589	 1.771	
## 
## sigma = 	1.000	1.795	1.000	1.777	2.374	1.000	1.000	1.396	2.415	1.000	
## 
## noiseSD = 	1	
```

```r
sim2 = ashsim(c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), rep(0.5, 7), rep(1/7, 7), niter = 20, 
    nsamp = 1000)
```

```
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8313	0.0821	0.0866	
## 
## mu = 		 0.000	-2.289	 2.243	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9976	0.0012	0.0012	
## 
## mu = 		0.0394	1.3313	1.5440	
## 
## sigma = 	1.485	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8143	0.0121	0.0296	0.0022	0.0558	0.0182	0.0162	0.0343	0.0123	0.0051	
## 
## mu = 		 0.000	-2.223	-2.223	-2.147	 2.167	-2.223	-2.223	 2.167	-2.223	 2.163	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## null?		 TRUE	FALSE	FALSE	FALSE	FALSE	 TRUE	FALSE	FALSE	FALSE	FALSE	
## 
## pi =		0.5957	0.0626	0.0408	0.0350	0.0424	0.0040	0.0224	0.1586	0.0288	0.0098	
## 
## mu = 		 0.0000	 1.6510	 1.6406	 1.6363	 1.6416	 0.0418	-1.8157	-1.8313	 1.6303	 1.5530	
## 
## sigma = 	1.000	1.002	1.001	1.001	1.001	1.665	1.000	1.000	1.001	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8405	0.0780	0.0815	
## 
## mu = 		 0.000	 2.225	-2.213	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9957	0.0015	0.0028	
## 
## mu = 		-0.0084	 0.9927	-2.0314	
## 
## sigma = 	1.443	1.227	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8283	0.0096	0.0097	0.0152	0.0016	0.0083	0.0175	0.0330	0.0292	0.0477	
## 
## mu = 		 0.0000	 2.1754	-2.1606	 2.1763	-0.2908	 2.1748	 2.1764	 2.1769	-2.1652	-2.1661	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.809	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6207	0.1657	0.0033	0.0427	0.0282	0.0393	0.0238	0.0037	0.0031	0.0694	
## 
## mu = 		 0.0000	-1.7095	-0.5418	 1.7240	 1.7233	 1.7239	-1.7026	 1.2924	-0.3061	 1.7245	
## 
## sigma = 	1.000	1.000	1.378	1.000	1.000	1.000	1.000	1.000	1.425	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8352	0.0918	0.0730	
## 
## mu = 		 0.000	-2.068	 2.257	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.9921	0.0053	0.0026	
## 
## mu = 		-0.0528	 2.2948	-0.0776	
## 
## sigma = 	1.440	1.000	1.437	
## 
## noiseSD = 	1	
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
## pi =		0.8165	0.0797	0.0172	0.0171	0.0171	0.0101	0.0083	0.0161	0.0025	0.0156	
## 
## mu = 		 0.000	 2.191	-1.992	-1.992	-1.992	-1.992	-1.992	-1.992	-1.925	-1.992	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## null?		 TRUE	 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	
## 
## pi =		0.5690	0.0026	0.0109	0.0774	0.0820	0.0395	0.0305	0.0672	0.1032	0.0177	
## 
## mu = 		 0.0000	 0.0706	-1.5442	 1.6736	-1.5497	-1.5489	 1.6642	 1.6723	-1.5501	 1.6559	
## 
## sigma = 	1.000	1.432	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7997	0.0963	0.1040	
## 
## mu = 		 0.000	-2.120	 2.092	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9958	0.0014	0.0028	
## 
## mu = 		 0.0331	 0.3717	-2.0026	
## 
## sigma = 	1.481	1.424	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7973	0.0231	0.0397	0.0433	0.0013	0.0114	0.0594	0.0094	0.0127	0.0024	
## 
## mu = 		 0.0000	-2.1360	 2.0689	 2.0689	-0.3628	 2.0689	-2.1360	-2.1360	 2.0689	-2.0599	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.811	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5570	0.0813	0.0807	0.0229	0.0045	0.0794	0.0521	0.0170	0.0096	0.0955	
## 
## mu = 		 0.000	 1.639	-1.676	-1.668	 1.451	-1.680	 1.638	-1.666	-1.652	 1.639	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8017	0.1021	0.0961	
## 
## mu = 		 0.000	-2.052	 2.244	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9935	0.0027	0.0038	
## 
## mu = 		-0.0223	 0.6729	-1.4795	
## 
## sigma = 	1.490	1.542	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8001	0.0292	0.0065	0.0014	0.0744	0.0345	0.0028	0.0082	0.0178	0.0252	
## 
## mu = 		 0.0000	-2.0242	-2.0234	 0.7454	 2.2762	-2.0242	-1.9630	-2.0239	 2.2762	-2.0242	
## 
## sigma = 	1.000	1.000	1.000	1.704	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## null?		 TRUE	FALSE	FALSE	 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	
## 
## pi =		0.5767	0.0089	0.0097	0.0024	0.0538	0.0115	0.2247	0.0269	0.0269	0.0585	
## 
## mu = 		 0.0000	 1.8460	-1.5457	-0.0914	 1.8655	 1.8515	-1.5783	 1.8605	 1.8605	 1.8660	
## 
## sigma = 	1.000	1.000	1.000	1.414	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8017	0.0982	0.1002	
## 
## mu = 		 0.000	-2.121	 2.206	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	 TRUE	
## 
## pi =		0.9970	0.0018	0.0012	
## 
## mu = 		 0.0163	-1.7430	 0.1715	
## 
## sigma = 	1.508	1.000	1.523	
## 
## noiseSD = 	1	
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
## pi =		0.8017	0.0961	0.0124	0.0022	0.0155	0.0203	0.0037	0.0141	0.0202	0.0138	
## 
## mu = 		 0.000	-2.135	 2.182	 2.740	 2.182	 2.182	 2.195	 2.182	 2.182	 2.182	
## 
## sigma = 	1.000	1.000	1.000	1.415	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5560	0.0049	0.0128	0.0393	0.0176	0.0242	0.0206	0.1633	0.0052	0.1560	
## 
## mu = 		 0.000	-1.438	 1.674	-1.640	 1.683	 1.688	-1.631	 1.706	 1.194	-1.657	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	2.629	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8024	0.0948	0.1028	
## 
## mu = 		 0.000	-2.027	 2.142	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9972	0.0013	0.0015	
## 
## mu = 		 0.0232	-1.0855	-1.3540	
## 
## sigma = 	1.471	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8035	0.0169	0.0288	0.0189	0.0201	0.0524	0.0307	0.0133	0.0130	0.0024	
## 
## mu = 		 0.000	 2.135	-2.046	 2.135	-2.046	 2.135	-2.046	 2.135	-2.045	 2.076	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.5520	0.0334	0.1196	0.0862	0.0034	0.0031	0.0278	0.0990	0.0107	0.0650	
## 
## mu = 		 0.0000	 1.7120	 1.7134	-1.5423	-0.7793	-0.3811	-1.5378	-1.5426	-1.5134	 1.7129	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.201	1.370	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8454	0.0722	0.0824	
## 
## mu = 		 0.000	 2.395	-2.217	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9974	0.0010	0.0016	
## 
## mu = 		-0.0053	-0.5007	 3.7952	
## 
## sigma = 	1.464	1.306	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8294	0.0059	0.0312	0.0341	0.0240	0.0024	0.0313	0.0065	0.0064	0.0287	
## 
## mu = 		 0.000	-2.156	-2.156	 2.310	-2.156	 2.202	 2.309	 2.298	 2.298	-2.156	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.5810	0.0918	0.0736	0.0197	0.0033	0.0270	0.0031	0.1004	0.0861	0.0140	
## 
## mu = 		 0.0000	-1.6610	-1.6605	 1.4322	 1.1288	-1.6592	 0.9324	 1.7255	 1.5042	-1.6574	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.730	1.000	1.740	1.174	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8271	0.0846	0.0883	
## 
## mu = 		 0.000	-2.234	 2.348	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9972	0.0013	0.0015	
## 
## mu = 		-0.008	-1.358	 2.365	
## 
## sigma = 	1.492	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8118	0.0280	0.0041	0.0189	0.0144	0.0092	0.0018	0.0197	0.0425	0.0497	
## 
## mu = 		 0.0000	-2.1559	-2.0704	-2.1535	-2.1518	-2.1462	-0.2383	-2.1537	 2.3102	 2.3103	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.866	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6229	0.0199	0.1038	0.0657	0.0028	0.0646	0.0650	0.0280	0.0210	0.0063	
## 
## mu = 		 0.000	-1.692	 1.927	-1.718	 1.463	 1.926	-1.718	-1.700	-1.693	 1.888	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7954	0.1009	0.1037	
## 
## mu = 		 0.000	-2.193	 2.061	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
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
## mu = 		 0.0058	-2.0993	-1.6427	
## 
## sigma = 	1.506	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7906	0.0293	0.0037	0.0174	0.0540	0.0014	0.0063	0.0026	0.0714	0.0233	
## 
## mu = 		 0.0000	 2.0376	-2.1651	-2.1922	 2.0376	-0.4732	-2.1894	-2.0949	-2.1925	 2.0376	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.864	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5253	0.0301	0.0030	0.0029	0.1195	0.0194	0.0763	0.0454	0.0967	0.0815	
## 
## mu = 		 0.0000	 1.5661	 0.6019	-0.4298	 1.5692	 1.5641	 1.5709	-1.7075	-1.7189	-1.7159	
## 
## sigma = 	1.000	1.000	1.337	1.562	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8030	0.0925	0.1046	
## 
## mu = 		 0.000	 2.185	-2.187	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9959	0.0025	0.0016	
## 
## mu = 		-0.0661	 1.9844	 1.5908	
## 
## sigma = 	1.506	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7957	0.0071	0.0231	0.0237	0.0259	0.0174	0.0604	0.0136	0.0026	0.0306	
## 
## mu = 		 0.000	 2.204	-2.131	 2.205	-2.131	-2.130	 2.205	-2.130	-2.017	-2.131	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.5665	0.0726	0.0817	0.0074	0.0470	0.0543	0.0769	0.0095	0.0759	0.0082	
## 
## mu = 		 0.000	 1.749	-1.700	-1.640	 1.747	-1.696	 1.749	-1.659	-1.699	-1.649	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8340	0.0862	0.0798	
## 
## mu = 		 0.000	 2.234	-2.246	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9963	0.0015	0.0022	
## 
## mu = 		0.0164	0.9566	1.8403	
## 
## sigma = 	1.478	1.197	1.000	
## 
## noiseSD = 	1	
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
## null?		 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	 TRUE	
## 
## pi =		0.8223	0.0103	0.0085	0.0520	0.0111	0.0045	0.0189	0.0048	0.0659	0.0018	
## 
## mu = 		 0.0000	 2.1750	-2.1941	-2.2139	 2.1761	-2.1392	-2.2070	 2.1291	 2.1881	-0.0817	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.916	
## 
## noiseSD = 	1	
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
## pi =		0.5939	0.1052	0.0211	0.0745	0.0050	0.1053	0.0030	0.0050	0.0292	0.0579	
## 
## mu = 		 0.0000	-1.7285	-1.7088	 1.6922	-1.5495	 1.6956	-0.2605	-1.5496	 1.6809	-1.7220	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.536	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8380	0.0685	0.0935	
## 
## mu = 		 0.000	-2.491	 2.243	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9912	0.0055	0.0033	
## 
## mu = 		 0.0626	-3.0217	 2.1740	
## 
## sigma = 	1.465	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8246	0.0321	0.0105	0.0175	0.0099	0.0235	0.0541	0.0150	0.0103	0.0027	
## 
## mu = 		 0.000	 2.191	-2.414	 2.190	-2.413	 2.190	-2.424	 2.190	 2.190	 2.093	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## null?		 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	 TRUE	FALSE	FALSE	
## 
## pi =		0.6689	0.0682	0.0565	0.0054	0.0030	0.0418	0.1177	0.0028	0.0202	0.0155	
## 
## mu = 		 0.0000	-2.0380	-2.0337	 1.6780	 0.5325	 1.8424	 1.8464	-0.0398	 1.8357	-1.9879	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.398	1.000	1.000	1.520	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8268	0.0816	0.0916	
## 
## mu = 		 0.000	-2.296	 2.278	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9954	0.0030	0.0016	
## 
## mu = 		0.0183	2.2762	1.3269	
## 
## sigma = 	1.495	1.000	1.103	
## 
## noiseSD = 	1	
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
## pi =		0.8096	0.0342	0.0046	0.0401	0.0277	0.0181	0.0194	0.0020	0.0036	0.0407	
## 
## mu = 		 0.000	-2.185	-2.186	 2.217	-2.185	 2.217	-2.185	-3.091	-2.195	 2.217	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.709	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5703	0.0370	0.0866	0.0784	0.0386	0.0036	0.1312	0.0329	0.0165	0.0049	
## 
## mu = 		 0.0000	 1.7133	 1.7215	 1.7205	-1.6606	-0.4298	-1.6643	-1.6595	 1.7011	-1.9981	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.617	1.000	1.000	1.000	2.229	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8322	0.1001	0.0677	
## 
## mu = 		 0.000	-2.141	 2.233	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9973	0.0015	0.0012	
## 
## mu = 		-0.0718	-1.8016	-1.4240	
## 
## sigma = 	1.46	1.00	1.00	
## 
## noiseSD = 	1	
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
## pi =		0.8225	0.0078	0.0248	0.0045	0.0067	0.0601	0.0224	0.0051	0.0310	0.0151	
## 
## mu = 		 0.000	 2.181	-2.108	 2.133	 2.175	-2.111	 2.196	-2.077	 2.198	-2.106	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.5617	0.0700	0.0612	0.0416	0.1784	0.0577	0.0091	0.0133	0.0027	0.0042	
## 
## mu = 		 0.000	-1.633	-1.633	-1.631	 1.611	-1.632	-1.622	 1.577	-0.852	 1.153	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.207	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8211	0.0689	0.1100	
## 
## mu = 		 0.000	 2.373	-2.346	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9979	0.0010	0.0011	
## 
## mu = 		-0.0953	 1.5488	 1.0711	
## 
## sigma = 	1.531	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8140	0.0155	0.0684	0.0301	0.0017	0.0206	0.0212	0.0059	0.0130	0.0095	
## 
## mu = 		 0.0000	 2.3302	-2.3359	 2.3435	-0.5335	 2.3361	-2.3202	 2.2897	-2.3124	-2.3045	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.961	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6193	0.0094	0.0034	0.0030	0.1658	0.0554	0.0037	0.0996	0.0232	0.0173	
## 
## mu = 		 0.0000	-1.9151	 0.6756	-1.2142	 1.7819	-1.9544	 0.9077	-1.9773	-1.9407	-1.9347	
## 
## sigma = 	1.000	1.000	1.485	1.484	1.002	1.000	1.356	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7848	0.1056	0.1097	
## 
## mu = 		 0.000	-2.158	 2.203	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9957	0.0023	0.0019	
## 
## mu = 		-0.0021	-1.7554	-1.5348	
## 
## sigma = 	1.541	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7735	0.0235	0.0023	0.0281	0.0115	0.0080	0.0131	0.0459	0.0549	0.0392	
## 
## mu = 		 0.000	-2.111	-1.994	-2.111	-2.111	-2.111	 2.173	 2.173	 2.173	-2.111	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.5194	0.0183	0.0210	0.0822	0.1888	0.0205	0.0371	0.0427	0.0069	0.0631	
## 
## mu = 		 0.000	-1.592	 1.771	 1.774	-1.635	 1.771	-1.616	 1.773	-1.438	 1.774	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8036	0.0986	0.0978	
## 
## mu = 		 0.000	-2.147	 2.137	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9971	0.0013	0.0016	
## 
## mu = 		 0.0117	-1.0558	-1.7760	
## 
## sigma = 	1.487	1.222	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8050	0.0176	0.0162	0.0406	0.0193	0.0178	0.0129	0.0063	0.0621	0.0022	
## 
## mu = 		 0.000	 2.144	-2.157	 2.144	 2.144	 2.144	-2.157	-2.157	-2.157	 2.032	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## null?		 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	 TRUE	FALSE	FALSE	FALSE	
## 
## pi =		0.5746	0.0063	0.0175	0.0506	0.0350	0.0707	0.0061	0.0099	0.0369	0.1923	
## 
## mu = 		 0.0000	-0.2824	 1.6644	 1.6898	 1.6831	 1.6964	-0.1756	-1.4754	 1.6841	-1.7384	
## 
## sigma = 	1.000	1.407	1.000	1.000	1.000	1.000	1.439	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8515	0.0752	0.0733	
## 
## mu = 		 0.000	 2.352	-2.219	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9969	0.0018	0.0013	
## 
## mu = 		 0.0432	-1.9302	 0.9767	
## 
## sigma = 	1.459	1.000	1.660	
## 
## noiseSD = 	1	
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
## pi =		0.8348	0.0095	0.0062	0.0413	0.0388	0.0173	0.0170	0.0155	0.0041	0.0155	
## 
## mu = 		 0.000	-2.153	-2.147	 2.269	 2.256	-2.154	-2.154	-2.155	 2.151	-2.154	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## null?		 TRUE	FALSE	FALSE	 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	
## 
## pi =		0.5721	0.0088	0.0962	0.0035	0.0944	0.0527	0.0355	0.0733	0.0162	0.0473	
## 
## mu = 		 0.0000	 1.3219	-1.6388	-0.0097	 1.7357	-1.6378	 1.4083	 1.4515	 1.3707	-1.6376	
## 
## sigma = 	1.000	1.000	1.000	1.544	1.196	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.8002	0.1057	0.0941	
## 
## mu = 		 0.000	-2.171	 2.085	
## 
## sigma = 	1	1	1	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9941	0.0033	0.0026	
## 
## mu = 		-0.0641	-2.0596	 1.8384	
## 
## sigma = 	1.488	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.8041	0.0061	0.0598	0.0117	0.0044	0.0040	0.0169	0.0511	0.0299	0.0120	
## 
## mu = 		 0.000	-2.137	 2.147	 2.140	-2.104	 2.059	-2.152	-2.155	-2.154	 2.141	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.5633	0.1828	0.0527	0.0202	0.0029	0.0095	0.0078	0.1291	0.0072	0.0245	
## 
## mu = 		 0.0000	-1.6867	-1.6802	 1.6801	 0.3023	 1.6655	 1.6533	 1.6882	 1.6464	 1.6814	
## 
## sigma = 	1.00	1.00	1.00	1.00	1.53	1.00	1.00	1.00	1.00	1.00	
## 
## noiseSD = 	1	
```

```r
sim3 = ashsim(c(-2, -1, 0, 1), c(2, 1.5, 1, 1), c(1/4, 1/4, 1/3, 1/6), niter = 20, 
    nsamp = 1000)
```

```
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7166	0.0830	0.2004	
## 
## mu = 		 0.000	 2.470	-3.055	
## 
## sigma = 	1.000	1.000	1.524	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9819	0.0146	0.0035	
## 
## mu = 		-0.3965	-6.1042	-4.8789	
## 
## sigma = 	1.881	1.272	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7052	0.0086	0.0024	0.0613	0.0313	0.0211	0.0145	0.0796	0.0201	0.0558	
## 
## mu = 		 0.000	-2.526	 2.298	-4.218	 2.413	 2.411	 2.410	-2.600	 2.411	-2.579	
## 
## sigma = 	1.000	1.000	1.000	1.717	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.3457	0.0643	0.0096	0.1227	0.0483	0.2848	0.0473	0.0045	0.0107	0.0622	
## 
## mu = 		 0.000	 2.091	-6.878	-1.391	 1.529	-2.171	 1.528	 1.451	-3.771	 1.546	
## 
## sigma = 	1.000	1.220	1.000	1.000	1.000	1.661	1.000	1.000	1.400	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7074	0.0846	0.2080	
## 
## mu = 		 0.000	 2.524	-3.113	
## 
## sigma = 	1.000	1.000	1.282	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9823	0.0164	0.0013	
## 
## mu = 		-0.4494	-5.2552	 1.5725	
## 
## sigma = 	1.92	1.00	1.00	
## 
## noiseSD = 	1	
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
## pi =		0.7016	0.0053	0.1222	0.0195	0.0301	0.0073	0.0050	0.0391	0.0067	0.0632	
## 
## mu = 		 0.000	-4.277	-2.961	-5.404	 2.491	 2.491	 2.491	 2.491	 2.491	-2.647	
## 
## sigma = 	1.000	1.267	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5498	0.0018	0.0406	0.0036	0.0695	0.0613	0.0477	0.1272	0.0569	0.0416	
## 
## mu = 		 0.0000	-0.4047	-4.9590	-1.9283	-2.4525	-2.4261	 2.1931	-2.5814	 2.1799	 2.1735	
## 
## sigma = 	1.000	1.806	1.039	1.134	1.034	1.030	1.000	1.057	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.6952	0.2279	0.0768	
## 
## mu = 		 0.000	-3.239	 2.453	
## 
## sigma = 	1.00	1.45	1.00	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9694	0.0284	0.0022	
## 
## mu = 		-0.5316	-5.3386	 1.3777	
## 
## sigma = 	1.945	1.028	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6850	0.0890	0.0286	0.0129	0.0154	0.0801	0.0051	0.0519	0.0217	0.0102	
## 
## mu = 		 0.000	-4.244	 2.438	 2.440	 2.438	-2.508	-3.833	-2.505	 2.438	-2.955	
## 
## sigma = 	1.000	1.386	1.000	1.000	1.000	1.000	1.485	1.000	1.000	1.051	
## 
## noiseSD = 	1	
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
## pi =		0.4427	0.0468	0.0151	0.2176	0.0517	0.0863	0.0074	0.0226	0.0660	0.0438	
## 
## mu = 		 0.000	 1.895	-1.733	-3.289	 1.915	-1.808	 1.736	 1.818	-1.786	 1.884	
## 
## sigma = 	1.000	1.049	1.000	1.681	1.069	1.000	1.000	1.000	1.000	1.038	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7256	0.2063	0.0681	
## 
## mu = 		 0.000	-3.059	 2.761	
## 
## sigma = 	1.000	1.432	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9768	0.0202	0.0030	
## 
## mu = 		-0.4273	-5.1748	-6.5295	
## 
## sigma = 	1.867	1.000	1.708	
## 
## noiseSD = 	1	
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
## pi =		0.7164	0.0107	0.0042	0.0297	0.0037	0.0451	0.0343	0.1251	0.0233	0.0075	
## 
## mu = 		 0.000	 2.720	-5.117	 2.722	-4.808	-2.550	-4.909	-2.621	 2.722	 2.717	
## 
## sigma = 	1.000	1.000	2.056	1.000	2.057	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5215	0.0652	0.1218	0.0073	0.0605	0.0194	0.1414	0.0029	0.0552	0.0047	
## 
## mu = 		 0.000	 2.204	-2.004	 1.831	 2.183	-1.920	-3.416	 1.516	-1.961	 1.691	
## 
## sigma = 	1.000	1.134	1.000	1.000	1.127	1.000	1.697	1.004	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.6877	0.0774	0.2349	
## 
## mu = 		 0.000	 2.776	-2.983	
## 
## sigma = 	1.000	1.000	1.475	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9911	0.0080	0.0009	
## 
## mu = 		-0.5587	-6.9787	 7.3635	
## 
## sigma = 	1.995	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6835	0.0038	0.0135	0.0370	0.0420	0.0393	0.0753	0.0844	0.0204	0.0010	
## 
## mu = 		 0.000	-4.022	-2.804	-4.858	 2.656	 2.656	-2.674	-2.695	-2.653	 7.304	
## 
## sigma = 	1.000	1.616	1.000	1.682	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
##  [1] 0.577016 0.010440 0.046455 0.026767 0.009882 0.207949 0.000000
##  [8] 0.032880 0.059247 0.029365
##  [1] -0.05416  1.26180 -2.49394  1.08537 -1.69689 -2.42759  0.12370
##  [8]  2.81437 -2.60651  2.58804
##  [1] 1.3020 0.5000 0.9910 2.6236 0.5774 2.1848 1.2240 0.7605 1.0055 0.8760
##           [,1]     [,2] [,3]
##  [1,] 0.990000  0.00000  0.5
##  [2,] 0.001111  0.67400  0.5
##  [3,] 0.001111 -1.17135  0.5
##  [4,] 0.001111  3.01095  0.5
##  [5,] 0.001111 -0.78096  0.5
##  [6,] 0.001111 -2.80793  0.5
##  [7,] 0.001111  0.02245  0.5
##  [8,] 0.001111  1.74920  0.5
##  [9,] 0.001111 -1.23313  0.5
## [10,] 0.001111  1.72941  0.5
##       [,1] [,2]
##  [1,] -Inf  Inf
##  [2,] -Inf  Inf
##  [3,] -Inf  Inf
##  [4,] -Inf  Inf
##  [5,] -Inf  Inf
##  [6,] -Inf  Inf
##  [7,] -Inf  Inf
##  [8,] -Inf  Inf
##  [9,] -Inf  Inf
## [10,] -Inf  Inf
##       [,1] [,2]
##  [1,]  0.5  Inf
##  [2,]  0.5  Inf
##  [3,]  0.5  Inf
##  [4,]  0.5  Inf
##  [5,]  0.5  Inf
##  [6,]  0.5  Inf
##  [7,]  0.5  Inf
##  [8,]  0.5  Inf
##  [9,]  0.5  Inf
## [10,]  0.5  Inf
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7058	0.0832	0.2110	
## 
## mu = 		 0.000	 2.479	-2.982	
## 
## sigma = 	1.000	1.000	1.381	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9780	0.0195	0.0025	
## 
## mu = 		-0.4089	-5.3505	-5.2495	
## 
## sigma = 	1.857	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6926	0.0053	0.0041	0.0194	0.0145	0.0344	0.0349	0.1605	0.0278	0.0066	
## 
## mu = 		 0.000	-2.776	-2.760	 2.432	-2.543	 2.432	-5.146	-2.543	 2.432	 2.432	
## 
## sigma = 	1	1	1	1	1	1	1	1	1	1	
## 
## noiseSD = 	1	
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
## pi =		0.4258	0.1112	0.0743	0.0654	0.1614	0.0050	0.0054	0.0750	0.0058	0.0708	
## 
## mu = 		 0.000	 1.839	-1.795	 1.828	-1.809	 2.066	 1.741	-3.982	 2.487	-2.422	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.807	1.000	1.429	1.682	1.745	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7392	0.1986	0.0621	
## 
## mu = 		 0.000	-3.064	 2.548	
## 
## sigma = 	1.000	1.213	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9829	0.0018	0.0153	
## 
## mu = 		-0.4866	 0.5565	-4.8441	
## 
## sigma = 	1.845	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7341	0.0059	0.0871	0.0054	0.0416	0.0126	0.0084	0.0116	0.0888	0.0045	
## 
## mu = 		 0.000	 2.535	-3.495	 2.528	 2.581	-2.530	 2.550	-2.528	-2.684	-3.590	
## 
## sigma = 	1.000	1.000	1.323	1.000	1.000	1.000	1.000	1.000	1.000	1.400	
## 
## noiseSD = 	1	
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
## pi =		0.5291	0.2083	0.0057	0.0032	0.0850	0.0146	0.0242	0.1000	0.0220	0.0081	
## 
## mu = 		 0.0000	-2.8721	-1.9375	-0.2805	 2.0504	-1.8495	 1.7612	-1.9768	 1.7510	 1.6025	
## 
## sigma = 	1.000	1.488	1.360	1.901	1.153	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7390	0.0669	0.1940	
## 
## mu = 		 0.000	 2.504	-3.225	
## 
## sigma = 	1.000	1.000	1.404	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9566	0.0013	0.0421	
## 
## mu = 		-0.310	 1.109	-4.917	
## 
## sigma = 	1.780	1.000	1.044	
## 
## noiseSD = 	1	
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
## pi =		0.7246	0.0593	0.0136	0.0045	0.0070	0.0322	0.0656	0.0112	0.0746	0.0074	
## 
## mu = 		 0.000	-2.540	 2.296	 3.590	 2.281	 2.342	-2.670	 2.291	-4.219	 2.282	
## 
## sigma = 	1.000	1.000	1.000	1.286	1.000	1.000	1.000	1.000	1.358	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4757	0.1367	0.0051	0.0520	0.0147	0.0889	0.0704	0.0866	0.0574	0.0124	
## 
## mu = 		 0.0000	-3.7568	 0.6514	 1.5256	 3.2724	-1.9387	 1.5493	-1.8864	 1.5317	-1.7811	
## 
## sigma = 	1.000	1.459	1.995	1.000	1.231	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7057	0.0789	0.2155	
## 
## mu = 		 0.000	 2.362	-3.036	
## 
## sigma = 	1.000	1.000	1.496	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9801	0.0015	0.0184	
## 
## mu = 		-0.4512	 1.1622	-5.8757	
## 
## sigma = 	1.884	1.000	1.237	
## 
## noiseSD = 	1	
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
## pi =		0.6954	0.0063	0.0070	0.0026	0.0656	0.0680	0.0266	0.0586	0.0089	0.0610	
## 
## mu = 		 0.000	 2.281	 2.283	 2.313	-2.626	-2.628	-2.586	-4.297	 2.287	 2.299	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.711	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4265	0.1521	0.1364	0.0171	0.0584	0.0120	0.0072	0.0200	0.1115	0.0588	
## 
## mu = 		 0.000	 1.712	-3.131	-6.222	 1.573	-1.711	-1.685	-1.740	-1.926	-1.823	
## 
## sigma = 	1.000	1.071	1.324	1.111	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7159	0.2147	0.0694	
## 
## mu = 		 0.000	-3.151	 2.518	
## 
## sigma = 	1.000	1.388	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9650	0.0012	0.0338	
## 
## mu = 		-0.4356	 1.2635	-5.0162	
## 
## sigma = 	1.852	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7025	0.0030	0.0807	0.0132	0.0617	0.0559	0.0195	0.0515	0.0041	0.0080	
## 
## mu = 		 0.000	 3.588	-4.123	 2.390	-2.557	-2.496	-2.478	 2.392	-3.831	 2.389	
## 
## sigma = 	1.000	1.494	1.340	1.000	1.000	1.000	1.000	1.000	1.506	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4334	0.0208	0.0101	0.0991	0.0076	0.0055	0.2620	0.0069	0.0333	0.1214	
## 
## mu = 		 0.000	 1.724	 1.690	 1.770	 2.918	 1.838	-2.727	-2.266	 1.735	-1.721	
## 
## sigma = 	1.000	1.000	1.000	1.002	1.547	1.886	1.755	2.123	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.6963	0.2185	0.0852	
## 
## mu = 		 0.000	-3.019	 2.513	
## 
## sigma = 	1.000	1.393	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9918	0.0068	0.0014	
## 
## mu = 		-0.5031	-6.7629	 1.2407	
## 
## sigma = 	1.977	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6909	0.0795	0.0158	0.0206	0.0036	0.0467	0.0167	0.0196	0.0891	0.0175	
## 
## mu = 		 0.000	-2.719	 2.466	 2.464	-4.087	-4.212	 2.463	 2.465	-2.719	 2.463	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.722	1.690	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4283	0.0163	0.1005	0.0046	0.1305	0.0751	0.0499	0.0030	0.1526	0.0393	
## 
## mu = 		 0.0000	-1.6733	-3.2478	 1.7536	-2.5296	-1.8908	-1.7922	-0.3966	 1.9210	 1.7099	
## 
## sigma = 	1.000	1.000	1.980	1.388	1.238	1.042	1.000	1.992	1.098	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7359	0.0602	0.2039	
## 
## mu = 		 0.000	 2.419	-2.941	
## 
## sigma = 	1.000	1.000	1.424	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9739	0.0246	0.0015	
## 
## mu = 		-0.4524	-5.1676	 1.6149	
## 
## sigma = 	1.746	1.206	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7236	0.0037	0.0463	0.0026	0.0120	0.0549	0.0084	0.0185	0.0643	0.0655	
## 
## mu = 		 0.000	 2.351	 2.370	 2.292	 2.369	-2.450	-2.402	-2.421	-4.046	-2.457	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.519	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5192	0.0102	0.0775	0.0401	0.0022	0.0492	0.0790	0.0774	0.1369	0.0082	
## 
## mu = 		 0.000	-1.699	-2.059	-1.818	 0.458	-1.852	-4.186	-2.057	 1.898	-1.681	
## 
## sigma = 	1.000	1.000	1.033	1.000	1.480	1.000	1.402	1.030	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.6768	0.0734	0.2499	
## 
## mu = 		 0.000	 2.538	-3.143	
## 
## sigma = 	1.000	1.000	1.427	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9811	0.0164	0.0025	
## 
## mu = 		-0.6316	-5.4708	-5.1786	
## 
## sigma = 	2.014	1.093	1.055	
## 
## noiseSD = 	1	
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
## pi =		0.6659	0.0034	0.0403	0.0060	0.0933	0.0367	0.0122	0.0978	0.0247	0.0197	
## 
## mu = 		 0.000	-3.686	 2.483	 2.481	-2.531	-2.521	 2.483	-4.080	-2.518	 2.483	
## 
## sigma = 	1.000	1.477	1.000	1.000	1.000	1.000	1.000	1.427	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.3804	0.0284	0.0032	0.0080	0.0116	0.0227	0.3182	0.1000	0.1168	0.0107	
## 
## mu = 		 0.000	 1.822	-1.098	 1.692	 1.742	-1.766	-2.667	-1.782	 1.913	 1.734	
## 
## sigma = 	1.000	1.034	1.994	1.000	1.000	1.000	1.815	1.000	1.081	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7152	0.2126	0.0722	
## 
## mu = 		 0.000	-3.033	 2.539	
## 
## sigma = 	1.000	1.543	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9722	0.0250	0.0028	
## 
## mu = 		-0.4334	-5.2819	-7.1883	
## 
## sigma = 	1.823	1.000	2.108	
## 
## noiseSD = 	1	
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
## pi =		0.7047	0.0188	0.0074	0.0307	0.0185	0.0502	0.0789	0.0040	0.0666	0.0202	
## 
## mu = 		 0.000	-2.418	 2.495	-2.420	 2.495	 2.495	-2.438	-6.148	-4.245	-2.419	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	2.499	1.383	1.000	
## 
## noiseSD = 	1	
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
## null?		 TRUE	FALSE	FALSE	FALSE	FALSE	FALSE	 TRUE	FALSE	FALSE	FALSE	
## 
## pi =		0.5219	0.0107	0.0024	0.0583	0.0544	0.0143	0.0017	0.1117	0.1090	0.1155	
## 
## mu = 		 0.0000	 2.0977	-1.3454	 2.1279	 2.1258	 2.1046	 0.0131	-2.0130	-1.9997	-3.8101	
## 
## sigma = 	1.000	1.000	1.863	1.000	1.000	1.000	1.753	1.000	1.000	1.746	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7124	0.1988	0.0888	
## 
## mu = 		 0.000	-3.165	 2.363	
## 
## sigma = 	1.000	1.233	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9793	0.0194	0.0012	
## 
## mu = 		-0.4252	-4.7899	 0.7291	
## 
## sigma = 	1.899	1.131	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7092	0.0141	0.0117	0.0695	0.0251	0.1083	0.0116	0.0289	0.0022	0.0196	
## 
## mu = 		 0.000	 2.362	-4.545	-3.322	 2.362	-2.773	-4.560	 2.362	 2.272	 2.362	
## 
## sigma = 	1.000	1.000	1.483	1.074	1.000	1.000	1.484	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.5215	0.0485	0.2405	0.0126	0.0050	0.0606	0.0327	0.0138	0.0583	0.0064	
## 
## mu = 		 0.000	 1.986	-2.767	-4.536	 1.919	 1.988	 1.983	 1.973	-2.082	 1.947	
## 
## sigma = 	1.000	1.000	1.376	1.810	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.6919	0.2222	0.0859	
## 
## mu = 		 0.000	-3.120	 2.479	
## 
## sigma = 	1.000	1.593	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9862	0.0124	0.0013	
## 
## mu = 		-0.5247	-6.4970	 1.4398	
## 
## sigma = 	1.978	2.032	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6963	0.0043	0.0115	0.0951	0.0763	0.0170	0.0081	0.0063	0.0119	0.0733	
## 
## mu = 		 0.000	-6.482	-2.582	-3.457	 2.455	-2.587	-6.936	-3.710	 2.455	-2.606	
## 
## sigma = 	1.000	3.146	1.000	1.121	1.000	1.000	1.000	1.121	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4805	0.1010	0.1249	0.0168	0.0440	0.0266	0.0557	0.0331	0.0832	0.0342	
## 
## mu = 		 0.000	-2.593	-3.124	-6.252	-1.692	 2.049	 2.056	-1.634	 2.059	-1.639	
## 
## sigma = 	1.000	1.363	1.323	2.154	1.000	1.000	1.000	1.000	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7128	0.1927	0.0945	
## 
## mu = 		 0.000	-3.136	 2.453	
## 
## sigma = 	1.000	1.531	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9715	0.0257	0.0028	
## 
## mu = 		-0.3504	-5.3751	 1.9166	
## 
## sigma = 	1.875	1.244	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.6997	0.0027	0.0117	0.0658	0.0085	0.0076	0.1163	0.0136	0.0563	0.0177	
## 
## mu = 		 0.000	-6.347	-5.186	 2.430	-4.901	-3.676	-2.359	 2.430	-3.701	 2.430	
## 
## sigma = 	1.000	2.759	1.000	1.000	1.000	1.317	1.000	1.000	1.317	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4799	0.0238	0.1130	0.0331	0.1232	0.0714	0.0010	0.0356	0.1043	0.0147	
## 
## mu = 		 0.000	 2.024	-1.854	 2.028	 2.042	-1.680	-9.727	-4.960	-3.082	-1.600	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.474	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.6982	0.2199	0.0819	
## 
## mu = 		 0.000	-2.958	 2.529	
## 
## sigma = 	1.000	1.367	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9883	0.0013	0.0105	
## 
## mu = 		-0.4892	 1.3742	-5.8012	
## 
## sigma = 	1.934	1.000	1.278	
## 
## noiseSD = 	1	
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
## pi =		0.6920	0.0985	0.0040	0.0361	0.0189	0.0796	0.0065	0.0221	0.0242	0.0183	
## 
## mu = 		 0.000	-2.595	-4.218	-4.377	 2.517	-2.594	-4.695	 2.517	 2.518	 2.517	
## 
## sigma = 	1.000	1.000	1.674	1.469	1.000	1.000	1.649	1.000	1.000	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4660	0.0072	0.1281	0.2085	0.0139	0.0666	0.0155	0.0645	0.0246	0.0052	
## 
## mu = 		 0.000	 1.711	-3.235	-2.125	 1.758	 1.822	-1.944	 2.163	 1.781	-1.837	
## 
## sigma = 	1.000	1.000	1.809	1.000	1.000	1.000	1.000	1.265	1.000	1.000	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7176	0.1931	0.0893	
## 
## mu = 		 0.000	-3.093	 2.512	
## 
## sigma = 	1.000	1.282	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9746	0.0242	0.0012	
## 
## mu = 		-0.3452	-4.7794	 0.7528	
## 
## sigma = 	1.882	1.000	1.006	
## 
## noiseSD = 	1	
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
## pi =		0.7057	0.0033	0.0107	0.0612	0.0842	0.0329	0.0769	0.0140	0.0059	0.0051	
## 
## mu = 		 0.000	-3.208	-2.440	-2.468	-3.714	-2.444	 2.497	 2.495	-3.112	-3.311	
## 
## sigma = 	1.000	1.315	1.000	1.000	1.310	1.000	1.000	1.000	1.160	1.246	
## 
## noiseSD = 	1	
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
## pi =		0.4751	0.0846	0.1094	0.0891	0.0441	0.0023	0.1414	0.0055	0.0451	0.0033	
## 
## mu = 		 0.000	-3.859	-1.712	 2.057	 1.916	 1.085	-2.274	 1.691	 1.919	-1.548	
## 
## sigma = 	1.000	1.306	1.000	1.101	1.021	1.395	1.293	1.000	1.024	1.420	
## 
## noiseSD = 	1	
## 
## 
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
## [1] "Warning: Posterior SDs not yet implemented for uniform components"
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
## pi =		0.7091	0.2278	0.0631	
## 
## mu = 		 0.000	-3.110	 2.469	
## 
## sigma = 	1.000	1.429	1.000	
## 
## noiseSD = 	1	
```

```
## Warning: Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization
## 
## Warning: Assuming known noise noiseSD =  1 . If needed rerun with noiseSD = NA to fit noiseSD.
## Warning: Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.
```

```
## Fitting preliminary model 
## Fitting final model 
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.9795	0.0182	0.0023	
## 
## mu = 		-0.5829	-5.5531	 1.1374	
## 
## sigma = 	1.906	1.397	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.7022	0.0035	0.0086	0.0037	0.0757	0.0054	0.0157	0.0115	0.1438	0.0299	
## 
## mu = 		 0.000	 2.341	 2.412	 2.349	-2.632	 2.388	 2.431	-2.513	-3.397	 2.463	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.562	1.000	
## 
## noiseSD = 	1	
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
## pi =		0.4205	0.0052	0.0078	0.0443	0.1519	0.1413	0.0859	0.0965	0.0310	0.0156	
## 
## mu = 		 0.000	-1.383	 2.738	 1.723	-2.195	-3.093	 1.735	-1.889	 1.717	-1.377	
## 
## sigma = 	1.000	1.000	1.519	1.000	1.336	1.962	1.000	1.220	1.000	1.000	
## 
## noiseSD = 	1	
```


Experiment with mixFdr to see how penalization effects fit. 

```r
betahat = sim1$betahat[[1]]
betahatsd = sim1$betahatsd[[1]]
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

res = data.frame(beta = x, mixfdr.P0 = t(cdf.mixfdr(fit.mixfdr.P0, x)), mixfdr = t(cdf.mixfdr(sim1$fit.mixfdr[[1]], 
    x)), ash.n = t(cdf.ash(sim1$fit.ash.n[[1]], x)$y), ash.fdr.n = t(cdf.ash(sim1$fit.ash.fdr.n[[1]], 
    x)$y))
```

```
## Error: no applicable method for 'mixcdf' applied to an object of class
## "normalmix"
```

```r

truth = t(cdf.ash(sim1$fit.ash.true[[1]], x)$y)

res.melt = melt(res, id.vars = c("beta"), variable.name = "Method")
```

```
## Error: could not find function "melt"
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




```r
plot_examples = function(sims) {
    sim1 = sims[[1]]
    sim2 = sims[[2]]
    sim3 = sims[[3]]
    len = length(sim1$beta[[1]])
    res = data.frame(beta = c(sim1$beta[[1]], sim2$beta[[1]], sim3$beta[[1]]), 
        Scenario = c(rep(1, len), rep(2, len), rep(3, len)))
    x = seq(-8, 6, length = 100)
    y = dnorm(x)
    dens1 = data.frame(x = x, dtrue = t(density(sim1$fit.ash.true[[1]], x)$y), 
        Scenario = 1)
    dens2 = data.frame(x = x, dtrue = t(density(sim2$fit.ash.true[[1]], x)$y), 
        Scenario = 2)
    dens3 = data.frame(x = x, dtrue = t(density(sim3$fit.ash.true[[1]], x)$y), 
        Scenario = 3)
    dens = rbind(dens1, dens2, dens3)
    ggplot(res) + facet_grid(. ~ Scenario) + geom_line(data = dens, aes(x, dtrue), 
        size = 1.2, alpha = 0.9, linetype = 1) + scale_x_continuous(limits = c(-6, 
        6))
}

pdf("figures/simABC_egdens.pdf", width = 6.5, height = 2)
plot_examples(list(sim1, sim2, sim3))
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
  rmse = function(x,y){sqrt(mean((x-y)^2))}
  get_rmse.ash = function(a,b){rmse(a$PosteriorMean,b)}
  get_rmse.mixfdr = function(a,b){rmse(a$effectSize,b)}
  plot_rmse = function(sims,inczero=FALSE,incbetahat=FALSE){
    res=list()
    
    for(i in 1:length(sims)){
      err.bayes = mapply(get_rmse.ash,sims[[i]]$fit.ash.true,sims[[i]]$beta)
      err.ash.n= mapply(get_rmse.ash,sims[[i]]$fit.ash.n,sims[[i]]$beta)
      err.ash.u = mapply(get_rmse.ash,sims[[i]]$fit.ash.u,sims[[i]]$beta)
      err.ash.hu = mapply(get_rmse.ash,sims[[i]]$fit.ash.hu,sims[[i]]$beta)
      err.mixfdr=mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr,sims[[i]]$beta)
      err.mixfdr.enull = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr,sims[[i]]$beta)
       err.mixfdr.J10 = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.J10,sims[[i]]$beta)
       #err.mixfdr.J100 = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.J100,sims[[i]]$beta)
      err.ash.fdr.n = mapply(get_rmse.ash,sims[[i]]$fit.ash.fdr.n,sims[[i]]$beta)
          
      err.betahat = mapply(rmse,sims[[i]]$betahat,sims[[i]]$beta)
      err.zero = unlist(lapply(sims[[i]]$beta,rmse,y=0))        
      res[[i]] = data.frame(Scenario=i,bayes=err.bayes,ash.normal=err.ash.n,ash.uniform=err.ash.u,ash.halfuniform=err.ash.hu,mixfdr=err.mixfdr,
                            mixfdr.enull = err.mixfdr.enull,
                            mixfdr.J10 = err.mixfdr.J10,
                            #mixfdr.J100 = err.mixfdr.J100,
                            ash.normal.fdr = err.ash.fdr.n)
      if(inczero){
        res[[i]]=data.frame(res[[i]],zero=err.zero)
      }
      if(incbetahat){
        res[[i]]=data.frame(res[[i]],betahat=err.betahat)
      }
    }
    require(reshape2)
    res.melt = melt(res, id.vars=c("Scenario","bayes"),variable.name="Method")    
    p=ggplot(data=res.melt,aes(bayes,value,colour=Method)) +geom_point(shape=16) +
         facet_grid(. ~ Scenario,scale="free_x") +
        geom_abline(colour = "black") +
        xlab("RMSE (Optimal Bayes Rule)") +
        ylab("RMSE (other method)")
    print(p +
          coord_equal(ratio=1))
  }
  pdf("figures/rmse_biplot.pdf",width=6.5,height=2)
    plot_rmse(list(sim1,sim2,sim3))
```

```
## Loading required package: reshape2
```

```r
  dev.off()
```

```
## pdf 
##   2
```

```r

#This figure not used in paper?
pdf("figures/rmse_biplot_withzerobetahat.pdf")
    plot_rmse(list(sim1,sim2,sim3),inczero=TRUE,incbetahat=TRUE)
  dev.off()
```

```
## pdf 
##   2
```





```r
get_rmse.ash = function(a, b) {
    rmse(a$PosteriorMean, b)
}
get_rmse.mixfdr = function(a, b) {
    rmse(a$effectSize, b)
}
plot_rmse_boxplot = function(sims, inczero = FALSE, incbetahat = FALSE, incmixfdr = FALSE) {
    res = list()
    
    for (i in 1:length(sims)) {
        err.bayes = mapply(get_rmse.ash, sims[[i]]$fit.ash.true, sims[[i]]$beta)
        err.ash.n = mapply(get_rmse.ash, sims[[i]]$fit.ash.n, sims[[i]]$beta)
        err.ash.u = mapply(get_rmse.ash, sims[[i]]$fit.ash.u, sims[[i]]$beta)
        err.ash.hu = mapply(get_rmse.ash, sims[[i]]$fit.ash.hu, sims[[i]]$beta)
        err.mixfdr = mapply(get_rmse.mixfdr, sims[[i]]$fit.mixfdr, sims[[i]]$beta)
        err.mixfdr.enull = mapply(get_rmse.mixfdr, sims[[i]]$fit.mixfdr, sims[[i]]$beta)
        err.mixfdr.J10 = mapply(get_rmse.mixfdr, sims[[i]]$fit.mixfdr.J10, sims[[i]]$beta)
        # err.mixfdr.J100 =
        # mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.J100,sims[[i]]$beta)
        err.ash.fdr.n = mapply(get_rmse.ash, sims[[i]]$fit.ash.fdr.n, sims[[i]]$beta)
        err.betahat = mapply(rmse, sims[[i]]$betahat, sims[[i]]$beta)
        err.zero = unlist(lapply(sims[[i]]$beta, rmse, y = 0))
        res[[i]] = data.frame(Scenario = i, ash.normal = err.ash.n/err.bayes, 
            ash.uniform = err.ash.u/err.bayes, ash.halfuniform = err.ash.hu/err.bayes, 
            ash.normal.fdr = err.ash.fdr.n/err.bayes)
        if (incmixfdr) {
            res[[i]] = data.frame(res[[i]], mixfdr = err.mixfdr/err.bayes, mixfdr.enull = err.mixfdr.enull/err.bayes, 
                mixfdr.J10 = err.mixfdr.J10/err.bayes)
            # mixfdr.J100 = err.mixfdr.J100/err.bayes)
        }
        if (inczero) {
            res[[i]] = data.frame(res[[i]], zero = err.zero/err.bayes)
        }
        if (incbetahat) {
            res[[i]] = data.frame(res[[i]], betahat = err.betahat/err.bayes)
        }
    }
    require(reshape2)
    res.melt = melt(res, id.vars = c("Scenario"), variable.name = "Method")
    ggplot(res.melt, aes(Method, value, color = Method)) + geom_boxplot() + 
        facet_grid(. ~ Scenario)
    
}
pdf("figures/rmse_boxplot.pdf", width = 6.5, height = 2)
plot_rmse_boxplot(list(sim1, sim2, sim3))
dev.off()
```

```
## pdf 
##   2
```

```r

pdf("figures/rmse_boxplot_extended.pdf", width = 6.5, height = 2)
plot_rmse_boxplot(list(sim1, sim2, sim3), TRUE, TRUE, TRUE)
dev.off()
```

```
## pdf 
##   2
```

```r

```



```r
plot_loglik_boxplot = function(sims) {
    res = list()
    for (i in 1:length(sims)) {
        loglik.bayes = mapply(get_loglik, sims[[i]]$fit.ash.true)
        loglik.ash.n = mapply(get_loglik, sims[[i]]$fit.ash.n)
        loglik.ash.u = mapply(get_loglik, sims[[i]]$fit.ash.u)
        loglik.ash.hu = mapply(get_loglik, sims[[i]]$fit.ash.hu)
        
        loglik.ash.fdr.n = mapply(get_loglik, sims[[i]]$fit.ash.fdr.n)
        res[[i]] = data.frame(Scenario = i, ash.normal = loglik.ash.n - loglik.bayes, 
            ash.uniform = loglik.ash.u - loglik.bayes, ash.halfuniform = loglik.ash.hu - 
                loglik.bayes, ash.normal.fdr = loglik.ash.fdr.n - loglik.bayes)
        
    }
    require(reshape2)
    res.melt = melt(res, id.vars = c("Scenario"), variable.name = "Method")
    ggplot(res.melt, aes(Method, value, color = Method)) + geom_boxplot() + 
        facet_grid(. ~ Scenario)
    
}
pdf("figures/loglik_boxplot.pdf", width = 6.5, height = 2)
plot_loglik_boxplot(list(sim1, sim2, sim3))
dev.off()
```

```
## pdf 
##   2
```

```r

```





```r
get_rmse.ash = function(a, b) {
    rmse(a$PosteriorMean, b)
}
get_rmse.mixfdr = function(a, b) {
    rmse(a$effectSize, b)
}
get_loglik.mixfdr = function(a, betahat, betahatsd) {
    loglik_conv(normalmix(a$pi, a$mu, a$sigma - 1), betahat, betahatsd)
}
cdf.mixfdr = function(a, x) {
    mixcdf(normalmix(a$pi, a$mu, a$sigma - 1), x)
}
plot_rmse_loglik_boxplot = function(sims) {
    res = list()
    
    for (i in 1:length(sims)) {
        err.bayes = mapply(get_rmse.ash, sims[[i]]$fit.ash.true, sims[[i]]$beta)
        err.ash.n = mapply(get_rmse.ash, sims[[i]]$fit.ash.n, sims[[i]]$beta)
        err.ash.u = mapply(get_rmse.ash, sims[[i]]$fit.ash.u, sims[[i]]$beta)
        err.ash.hu = mapply(get_rmse.ash, sims[[i]]$fit.ash.hu, sims[[i]]$beta)
        err.mixfdr = mapply(get_rmse.mixfdr, sims[[i]]$fit.mixfdr, sims[[i]]$beta)
        err.ash.fdr.n = mapply(get_rmse.ash, sims[[i]]$fit.ash.fdr.n, sims[[i]]$beta)
        
        res[[i]] = data.frame(Scenario = i, ash.normal = err.ash.n/err.bayes, 
            ash.uniform = err.ash.u/err.bayes, ash.halfuniform = err.ash.hu/err.bayes, 
            ash.normal.fdr = err.ash.fdr.n/err.bayes)
        
    }
    
    res2 = list()
    for (i in 1:length(sims)) {
        loglik.bayes = mapply(get_loglik, sims[[i]]$fit.ash.true)
        loglik.ash.n = mapply(get_loglik, sims[[i]]$fit.ash.n)
        loglik.ash.u = mapply(get_loglik, sims[[i]]$fit.ash.u)
        loglik.ash.hu = mapply(get_loglik, sims[[i]]$fit.ash.hu)
        loglik.mixfdr = mapply(get_loglik.mixfdr, sims[[i]]$fit.mixfdr, sims[[i]]$betahat, 
            sims[[i]]$betahatsd)
        loglik.ash.fdr.n = mapply(get_loglik, sims[[i]]$fit.ash.fdr.n)
        res2[[i]] = data.frame(Scenario = i, ash.normal = loglik.ash.n - loglik.bayes, 
            ash.uniform = loglik.ash.u - loglik.bayes, ash.halfuniform = loglik.ash.hu - 
                loglik.bayes, ash.normal.fdr = loglik.ash.fdr.n - loglik.bayes)
        
    }
    
    require(reshape2)
    res.melt = melt(res, id.vars = c("Scenario"), variable.name = "Method")
    res2.melt = melt(res2, id.vars = c("Scenario"), variable.name = "Method")
    res.melt$type = "RMSE (vs Bayes Rule)"
    res2.melt$type = "log(likelihood) (vs Bayes Rule)"
    ggplot(rbind(res.melt, res2.melt), aes(Method, value)) + geom_boxplot() + 
        facet_grid(type ~ Scenario, scale = "free_y")
    
}
pdf("figures/rmse_loglik_boxplot.pdf", width = 6.5, height = 5)
plot_rmse_loglik_boxplot(list(sim1, sim2, sim3))
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
plot(sim3$fit.ash.hu[[1]]$PosteriorMean, sim3$fit.mixfdr[[1]]$effectSize)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-101.png) 

```r
plot(sim3$fit.ash.fdr.n[[1]]$PosteriorMean, sim3$fit.mixfdr[[1]]$effectSize)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-102.png) 


problem seems to be that it overestimates pi0

```r
sim3$fit.mixfdr[[1]]$pi
```

```
## [1] 0.71661 0.08298 0.20041
```

[1] 0.71660654 0.08298064 0.20041282



```r
temp = mixFdr(sim3$betahat[[1]]/sim3$betahatsd[[1]], noiseSD = 1, theonull = TRUE, 
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

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) 

```
## 
## Fitted Model: J = 3 groups
## ----------------------------
## null?		 TRUE	FALSE	FALSE	
## 
## pi =		0.6159	0.1227	0.2614	
## 
## mu = 		 0.000	 2.200	-2.721	
## 
## sigma = 	1.000	1.000	1.639	
## 
## noiseSD = 	1	
```

```r
temp2 = mixFdr(sim3$betahat[[1]]/sim3$betahatsd[[1]], noiseSD = 1, theonull = TRUE, 
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
## pi =		0.7076	0.0279	0.0079	0.0279	0.0081	0.0044	0.0433	0.1178	0.0249	0.0303	
## 
## mu = 		 0.000	 2.413	-2.508	 2.413	 2.410	-3.718	-4.707	-2.665	 2.413	-2.530	
## 
## sigma = 	1.000	1.000	1.000	1.000	1.000	1.380	1.619	1.000	1.000	1.000	
## 
## noiseSD = 	1	
```

```r
plot(ecdf(sim3$beta[[1]]), xlim = c(-6, 6))
lines(x, cdf.mixfdr(temp, x), col = 3)
```

```
## Error: no applicable method for 'mixcdf' applied to an object of class
## "normalmix"
```

```r
lines(cdf.ash(sim3$fit.ash.n[[1]], x), col = 2)
lines(cdf.ash(sim3$fit.ash.hu[[1]], x), col = 2, lty = 2)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) 





```r
plot_LR = function(sims) {
    hist(unlist(lapply(sims$fit.ash.u, get_loglik)) - unlist(lapply(sims$fit.ash.n, 
        get_loglik)), xlab = "loglik difference", main = "loglik differences for nullbiased prior vs mle", 
        nclass = 10)
}

pdf("figures/logLR.pdf")
plot_LR(sim1)
plot_LR(sim2)
dev.off()
```

```
## pdf 
##   2
```


## Unused figures?


```r
plot_examples_withfit = function(sims) {
    sim1 = sims[[1]]
    sim2 = sims[[2]]
    sim3 = sims[[3]]
    len = length(sim1$beta[[1]])
    res = data.frame(beta = c(sim1$beta[[1]], sim2$beta[[1]], sim3$beta[[1]]), 
        Scenario = c(rep(1, len), rep(2, len), rep(3, len)))
    x = seq(-8, 6, length = 100)
    y = dnorm(x)
    dens1 = data.frame(x = x, dn = t(density(sim1$fit.ash.n[[1]], x)$y), du = t(density(sim1$fit.ash.u[[1]], 
        x)$y), dhu = t(density(sim1$fit.ash.hu[[1]], x)$y), Scenario = 1)
    dens2 = data.frame(x = x, dn = t(density(sim2$fit.ash.n[[1]], x)$y), du = t(density(sim2$fit.ash.u[[1]], 
        x)$y), dhu = t(density(sim2$fit.ash.hu[[1]], x)$y), Scenario = 2)
    dens3 = data.frame(x = x, dn = t(density(sim3$fit.ash.n[[1]], x)$y), du = t(density(sim3$fit.ash.u[[1]], 
        x)$y), dhu = t(density(sim3$fit.ash.hu[[1]], x)$y), Scenario = 3)
    dens = rbind(dens1, dens2, dens3)
    cbbPalette <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#000000", 
        "#D55E00", "#CC79A7")
    ggplot(res, aes(x = beta)) + facet_grid(. ~ Scenario) + geom_histogram(aes(y = ..density..), 
        binwidth = 0.5, colour = "black", fill = "white") + geom_line(data = dens, 
        aes(x, du, color = "Uniform"), size = 1.2, alpha = 0.9, linetype = 1) + 
        geom_line(data = dens, aes(x, dhu, color = "Half Uniform"), size = 1.2, 
            alpha = 0.9, linetype = 1) + geom_line(data = dens, aes(x, dn, color = "Normal"), 
        size = 1.2, alpha = 0.9, linetype = 1) + scale_colour_manual(name = "Method", 
        values = cbbPalette) + scale_x_continuous(limits = c(-4, 4))
}

pdf("figures/simABC_eg_withfit.pdf", width = 6.5, height = 2)
plot_examples_withfit(list(sim1, sim2, sim3))
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
cdf.mixfdr = function(a, x) {
    mixcdf(normalmix(a$pi, a$mu, a$sigma - 1), x)
}
plot_examples_cdf_withfit = function(sims) {
    sim1 = sims[[1]]
    sim2 = sims[[2]]
    sim3 = sims[[3]]
    len = length(sim1$beta[[1]])
    res = data.frame(beta = c(sim1$beta[[1]], sim2$beta[[1]], sim3$beta[[1]]), 
        Scenario = c(rep(1, len), rep(2, len), rep(3, len)))
    x = seq(-6, 6, length = 100)
    y = dnorm(x)
    cdf1 = data.frame(x = x, dn = t(cdf.ash(sim1$fit.ash.n[[1]], x)$y), du = t(cdf.ash(sim1$fit.ash.u[[1]], 
        x)$y), dhu = t(cdf.ash(sim1$fit.ash.hu[[1]], x)$y), dtrue = t(cdf.ash(sim1$fit.ash.true[[1]], 
        x)$y), dn.fdr = t(cdf.ash(sim1$fit.ash.fdr.n[[1]], x)$y), dm = t(cdf.mixfdr(sim1$fit.mixfdr[[1]], 
        x)), Scenario = 1)
    cdf2 = data.frame(x = x, dn = t(cdf.ash(sim2$fit.ash.n[[1]], x)$y), du = t(cdf.ash(sim2$fit.ash.u[[1]], 
        x)$y), dhu = t(cdf.ash(sim2$fit.ash.hu[[1]], x)$y), dtrue = t(cdf.ash(sim2$fit.ash.true[[1]], 
        x)$y), dn.fdr = t(cdf.ash(sim2$fit.ash.fdr.n[[1]], x)$y), dm = t(cdf.mixfdr(sim2$fit.mixfdr[[1]], 
        x)), Scenario = 2)
    cdf3 = data.frame(x = x, dn = t(cdf.ash(sim3$fit.ash.n[[1]], x)$y), du = t(cdf.ash(sim3$fit.ash.u[[1]], 
        x)$y), dhu = t(cdf.ash(sim3$fit.ash.hu[[1]], x)$y), dtrue = t(cdf.ash(sim3$fit.ash.true[[1]], 
        x)$y), dn.fdr = t(cdf.ash(sim3$fit.ash.fdr.n[[1]], x)$y), dm = t(cdf.mixfdr(sim3$fit.mixfdr[[1]], 
        x)), Scenario = 3)
    
    cdf = rbind(cdf1, cdf2, cdf3)
    cbbPalette <- c("#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#E69F00", 
        "#000000", "#D55E00", "#CC79A7")
    ggplot(res, aes(x = beta)) + facet_grid(. ~ Scenario) + geom_line(data = cdf, 
        aes(x, dtrue, color = "True"), size = 1.5, alpha = 0.9, linetype = 1) + 
        geom_line(data = cdf, aes(x, du, color = "Uniform"), size = 1, alpha = 0.9, 
            linetype = 2) + geom_line(data = cdf, aes(x, dhu, color = "Half Uniform"), 
        size = 1, alpha = 0.9, linetype = 2) + geom_line(data = cdf, aes(x, 
        dn, color = "Normal"), size = 1, alpha = 0.9, linetype = 2) + geom_line(data = cdf, 
        aes(x, dn.fdr, color = "Normal, null-biased"), size = 1, alpha = 0.9, 
        linetype = 2) + geom_line(data = cdf, aes(x, dm, color = "mixfdr"), 
        size = 1, alpha = 0.9, linetype = 2) + scale_colour_manual(name = "Method", 
        values = cbbPalette) + scale_x_continuous(limits = c(-6, 6))
}

pdf("figures/simABC_eg_cdf_withfit.pdf", width = 6.5, height = 2)
plot_examples_cdf_withfit(list(sim1, sim2, sim3))
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


