#the aim of this dsc is to assess the distribution of the logLR 
#value for ash under the null that all effects are 0

library(dscr)

dsc_znull=new.dsc("znull","dsc-znull-files")
source("add_scenarios.null.R")
source("add_methods.null.R")

source("score.R")
add_score(dsc_znull,score_logLR,"logLR","ash_output")


res=run_dsc(dsc_znull)
save(res,dsc_znull,file="dsc-znull-files/res.znull.RData")

