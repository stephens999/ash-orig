This repository contains ongoing work on Bayesian FDR and Adaptive Shrinkage.

# The ashr R Package 

The (in development) ashr package is in the `package` subdirectory.

There are two versions of the package, one which uses some C++ via RCpp, and one that does not.
The C++ version is faster, but even the regular R package is reasonably quick and
the extra speed is likely only required if using ash 
many times on large data sets. Since installing the C++ version requires more
dependencies (e.g. you need to have a C++ compiler installed etc) you might want to start with 
the non C++ version.

## Non C++ version

To install the pure-R version of the package, first download `package/ashr.no.cxx.tar.gz`.
You can download this either by cloning the whole repo, or by, for example, going to 
[https://github.com/stephens999/ash/blob/master/package/ashr.no.cxx.tar.gz](https://github.com/stephens999/ash/blob/master/package/ashr.no.cxx.tar.gz) and click "View Raw".
Your browser will likely download the file automatically and save it in some
default place (e.g. "Downloads"). 

After saving the file somewhere on your computer, you can install the package as follows: 
```
> install.packages("truncnorm")
> install.packages("/path/to/ashr.no.cxx.tar.gz",repos=NULL,type="source")
```
Here `/path/to/` indicates where on your computer you have saved the file.


## C++ version

We also have a faster version of the package, which exploits the very helpful Rcpp package
to interface R with C++.

To install the C++ version, first download `package/ashr.tar.gz`.
You can download this either by cloning the whole repo, or by, for example, going to
[https://github.com/stephens999/ash/blob/master/package/ashr.tar.gz](https://github.com/stephens999/ash/blob/master/package/ashr.tar.gz) and click "View Raw".
Your browser will likely download the file automatically and save it in some
default place (e.g. "Downloads"). 

After saving the file somewhere on your computer, you can install the package as follows:
```
> install.packages("truncnorm")
> install.packages("Rcpp")
> install.packages("/path/to/ashr.tar.gz",repos=NULL,type="source")
```
Here `/path/to/` indicates where on your computer you have saved the file.

## Running Adaptive SHrinkage


The main function in the ashr package is `ash`. To get minimal help:
```
> library("ashr")
> ?ash
```

To get more information download the repository, and look at `Rcode/intro.html` and `Rcode/readme.html`

(These html files are generated from `Rcode/intro.rmd` and `Rcode/readme.rmd` and can be
reproduced using the R package knitr. R studio provides a nice interface to accomplish this.
See for example http://jeromyanglim.blogspot.com/2012/05/getting-started-with-r-markdown-knitr.html)

If you are interested in contributing to this project, contact mstephens@uchicago.edu

