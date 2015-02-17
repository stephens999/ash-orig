library(dscr)

source("scenarios.R")
source("methods.R")
source("score.R")
res=run_dsc(scenarios,methods,score)
save(res,file="res.RData")



