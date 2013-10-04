#These commands were used to set up the ashr package

setwd("~/Documents/git/ash/package")
library(devtools)
create("~/Documents/git/ash/package/ashr")

Then I had to add 
export(ash,EMest)
to NAMESPACE file

#Then I put ash.R and mix.R in the R directory
#NOte: I had to remove the source("") and library("") commands from these files
#to get them to work

setwd("~/Documents/git/ash/package/ashr")
# not run build(binary=TRUE)
install("../ashr")
ash
EMest
works!


# had problems with these,
R CMD check ashr
R CMD INSTALL -l "/Applications/RStudio.app/Contents/Resources/R/library" ashr
install.packages("/Applications/RStudio.app/Contents/Resources/R/library/ashr")
- that lead to ash not actually appearing - not sure why!