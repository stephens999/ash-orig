#put all files (R and cpp files) in a directory and cd into that directory
#in R:
library(Rcpp)
Rcpp.package.skeleton("ashr", path=".", code_files=c("ash.oldversion.R", "ash.R", "mix.R"), cpp_files="MixEM.cpp", example_code = FALSE, attributes = TRUE)
#in the terminal
#cd ..
#tar -pczf ashr.tar.gz ashr
#in R
install.packages("ash/package/ashr.tar.gz",repos=NULL,type="source")
