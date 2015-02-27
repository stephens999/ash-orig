[![Build Status](https://travis-ci.org/daichaoxing/ash.svg?branch=travis)](https://travis-ci.org/daichaoxing/ash)

This repository contains ongoing work on Bayesian FDR and Adaptive Shrinkage.
Because the repo was getting rather out of hand I moved the package
to its own repo. So:
the ashr R package is now [here](https://github.com/stephens999/ashr)

To install the ashr package first you need to install devtools
```
install.packages("devtools")
library(devtools)
install_github("stephens999/ashr")
```

## Running Adaptive Shrinkage


The main function in the ashr package is `ash`. To get minimal help:
```
> library(ashr)
> ?ash
```

To get more information download the repository, and look at `Rcode/intro.html` and `Rcode/readme.html`

(These html files are generated from `Rcode/intro.rmd` and `Rcode/readme.rmd` and can be
reproduced using the R package knitr. R studio provides a nice interface to accomplish this.
See for example http://jeromyanglim.blogspot.com/2012/05/getting-started-with-r-markdown-knitr.html)

If you are interested in contributing to this project, contact mstephens@uchicago.edu

