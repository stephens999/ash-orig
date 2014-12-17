#install.packages("/Users/stephens/Dropbox/Documents/git/dscr_0.1.tar.gz",repos=NULL,type="source")

library("plyr")
library("ashr")
library("reshape2")
library(dscr)

source("datamaker.R")
source("method.R")
source("scenario.R")
source("score.R")

res=run_dsc(scenarios,methods,score)


head(res)

aggregate(td~method+scenario,res,mean)
aggregate(fd~method+scenario,res,mean)
res$fdr = res$fd/(res$fd+res$td)
aggregate(fdr~method+scenario,res,mean)
