Test turboEM
========================================================

First we have to define a fixpoint function and an objective function
for turboEM to work with

```r
normalize = function(x) {
    return(x/sum(x))
}

fixpoint = function(pi, matrix_lik, prior) {
    m = t(pi * t(matrix_lik))  # matrix_lik is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    classprob = m/m.rowsum  #an n by k matrix
    pinew = normalize(colSums(classprob) + prior - 1)
    return(pinew)
}

negpenloglik = function(pi, matrix_lik, prior) {
    m = t(pi * t(matrix_lik))  # matrix_lik is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    loglik = sum(log(m.rowsum))
    priordens = sum((prior - 1) * log(pi))
    return(-(loglik + priordens))
}

# Note this doesn't impose the constraint that they have to sum to 1...
# might be a problem?
pconstr <- function(par) {
    return(all(0 <= par & par <= 1))
}
```


Now simulate some test data

```r
set.seed(100)
sampsize = 10000
sd = c(1, 1.1, 1.2)
z = rnorm(sampsize, 0, sample(sd, sampsize, replace = TRUE))
lik = t(vapply(z, dnorm, sd, sd = sd))
prior = c(1, 1, 1)
pi.init = c(0.3, 0.2, 0.5)
```



```r
library(turboEM)
```

```
## Loading required package: numDeriv
## Loading required package: quantreg
## Loading required package: SparseM
## 
## Attaching package: 'SparseM'
## 
## The following object is masked from 'package:base':
## 
##     backsolve
## 
## Loading required package: foreach
## 
## Attaching package: 'turboEM'
## 
## The following objects are masked from 'package:numDeriv':
## 
##     grad, hessian
```

```r
res = turboem(par = pi.init, control.run = list(convtype = "objfn", tol = 1e-05), 
    fixptfn = fixpoint, objfn = negpenloglik, pconstr = pconstr, method = c("em", 
        "squarem", "pem", "decme", "qn"), matrix_lik = lik, prior = prior)
```

```
## Warning: NaNs produced
## Warning: NaNs produced
## Warning: NaNs produced
```

```r
options(digits = 13)
res
```

```
##    method    value.objfn  itr fpeval objfeval convergence elapsed.time
## 1      em 15161.96849170 1500   1500     1501       FALSE        4.118
## 2 squarem 15161.69882417   18     35       19        TRUE        0.119
## 3     pem 15161.70092168   13     36      102        TRUE        0.204
## 5      qn 15161.70340886    8     12       16        TRUE        0.044
## 
## Acceleration scheme 4 (decme) failed
```

```r
library(devtools)
devtools::load_all("../package/ashr")
```

```
## Loading ashr
## Loading required namespace: truncnorm
## Loading required namespace: SQUAREM
## Loading required namespace: Rcpp
## Loading required package: truncnorm
## Loading required package: SQUAREM
## Loading required package: Rcpp
```

```r
system.time(res2 <- mixEM(lik, prior, pi.init))
```

```
##    user  system elapsed 
##   0.080   0.003   0.085
```

```r
length(res2$B)
```

```
## [1] 1
```

```r
res2$B[length(res2$B)]
```

```
## [1] 15161.69882417
```


The number of iterations for the ash implementation is 1 and the objective achieved
is 1.5161698824165 &times; 10<sup>4</sup>.

Now try a bigger sample, just comparing EM and squareEM

```r
set.seed(100)
sampsize = 5e+05
sd = c(1, 1.1, 1.2)
z = rnorm(sampsize, 0, sample(sd, sampsize, replace = TRUE))
lik = t(vapply(z, dnorm, sd, sd = sd))
prior = c(1, 1, 1)
pi.init = c(0.3, 0.2, 0.5)
res = turboem(par = pi.init, control.run = list(convtype = "objfn", tol = 1e-05), 
    fixptfn = fixpoint, objfn = negpenloglik, pconstr = pconstr, method = c("em", 
        "squarem"), matrix_lik = lik, prior = prior)
options(digits = 13)
res
```

```
##    method    value.objfn  itr fpeval objfeval convergence elapsed.time
## 1      em 758092.1615915 1500   1500     1501       FALSE      292.048
## 2 squarem 758092.0206260   23     45       24        TRUE        8.711
```


Now run the squarem package directly

```r
library(SQUAREM)
system.time(res3 <- squarem(par = pi.init, fixptfn = fixpoint, objfn = negpenloglik, 
    matrix_lik = lik, prior = prior))
```

```
##    user  system elapsed 
##   7.494   0.949   8.459
```


It seems from this limited assessment that i) the squareEM approach is the most effective,
and ii) the SQUAREM package implementation is faster then the turboem implementation.

