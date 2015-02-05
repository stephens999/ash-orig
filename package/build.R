#put all files (R and cpp files) in a directory and cd into that directory
#in R:
library(Rcpp)
Rcpp.package.skeleton("ashr", path=".", code_files=c("ashutility.R", "ash.R", "mix.R"), cpp_files="MixEM.cpp", example_code = FALSE, attributes = TRUE)

#require(roxygen2)
#require(devtools)
#roxygenize("~/Documents/git/ash/package/ashr/R/ashr")
#document("~/Documents/git/ash/package/ashr/R/ashr")

#in the terminal
#tar -pczf ashr.tar.gz ashr
#sh create.ashr.no.cxx.sh
#in R
install.packages("~/Dropbox/Documents/git/ash/package/ashr.tar.gz",repos=NULL,type="source")
