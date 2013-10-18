#This file documents some of the learning steps I took as I set up the ashr package

#This used to create package skeleton
require(devtools)
create("~/Documents/git/ash/package/ashr")

#THEN: I had to add export(ash) to NAMESPACE file
#and  mv ash.R 
# mv mix.R 
#into the R directory

#NOte: I had to remove the source("") and library("") commands from the
#ash.R and mix.R files to get them to work

# process instructions; this als sets up NAMESPACE file
require(roxygen2)
roxygenize("~/Documents/git/ash/package/ashr")

#install from source: this uses R CMD INSTALL to install
# the ashr library in the usual directory (.Library)
#and loads it
require(devtools)
document("~/Documents/git/ash/package/ashr")
install("~/Documents/git/ash/package/ashr")
#I found I need to restart R after install to avoid errors in the help commands
library(ashr)
ash
?ash
#works!

#making a tar.gz version of ashr package
require(devtools)
build("~/Documents/git/ash/package/ashr",binary=FALSE)

#this installs the tar.gz package and puts it in the library
install.packages("package/ashr_0.1.tar.gz",repos=NULL,type="source")
library(ashr)
?ash
#works!

#Try installing from github; worked except for documentation - not sure why
require(devtools) 
install_github('ash',user='stephens999',subdir='package/ashr')
library(ashr)
ash
?ash


# had problems with these,
R CMD check ashr
R CMD INSTALL -l "/Applications/RStudio.app/Contents/Resources/R/library" ashr
install.packages("/Applications/RStudio.app/Contents/Resources/R/library/ashr")
- that lead to ash not actually appearing - not sure why!