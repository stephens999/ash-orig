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

penloglik = function(pi, matrix_lik, prior) {
    m = t(pi * t(matrix_lik))  # matrix_lik is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    loglik = sum(log(m.rowsum))
    priordens = sum((prior - 1) * log(pi))
    return(loglik + priordens)
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
    fixptfn = fixpoint, objfn = penloglik, pconstr = pconstr, method = c("em", 
        "squarem", "pem", "decme", "qn"), matrix_lik = lik, prior = prior)
options(digits = 13)
res
```

```
##    method     value.objfn  itr fpeval objfeval convergence elapsed.time
## 1      em -15161.96849170 1500   1500     1501       FALSE        4.634
## 2 squarem -15161.69884713   22     43       25        TRUE        0.163
## 3     pem -15161.90144876 1500   3010     3001       FALSE        9.081
## 5      qn -15161.70968009  913    917     1827        TRUE        6.255
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
## Loading required namespace: Rcpp
## Loading required package: truncnorm
## Loading required package: Rcpp
```

```r
system.time(res2 <- mixEM(lik, prior, pi.init))
```

```
##    user  system elapsed 
##   0.657   0.016   0.676
```

```r
length(res2$B)
```

```
## [1] 320
```

```r
res2$B[length(res2$B)]
```

```
## [1] -15162.041295
```


The number of iterations for the ash implementation is 320 and the objective achieved
is -1.5162041295002 &times; 10<sup>4</sup>.

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
    fixptfn = fixpoint, objfn = penloglik, pconstr = pconstr, method = c("em", 
        "squarem"), matrix_lik = lik, prior = prior)
options(digits = 13)
res
```

```
##    method     value.objfn  itr fpeval objfeval convergence elapsed.time
## 1      em -758092.1615915 1500   1500     1501       FALSE      309.421
## 2 squarem -758092.0206264  116    231      195        TRUE       49.537
```


