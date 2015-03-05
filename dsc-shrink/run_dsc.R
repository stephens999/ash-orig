library(dscr)

dsc_shrink=new.dsc("shrink","dsc-shrink-files")
source("addScenarios.R")
source("addMethods.R")

addOutputParser(dsc_shrink,"ash2beta",ash2beta_est,"ash_output","beta_est_output")
addOutputParser(dsc_shrink,"mixfdr2beta",mixfdr2beta_est,"mixfdr_output","beta_est_output")

addOutputParser(dsc_shrink,"ash2pi0",ash2pi0_est,"ash_output","pi0_est_output")
addOutputParser(dsc_shrink,"mixfdr2pi0",mixfdr2pi0_est,"mixfdr_output","pi0_est_output")


source("score.R")
addScore(dsc_shrink,score,"beta_err","beta_est_output")
addScore(dsc_shrink,score2,"pi0_score","pi0_est_output")


res=run_dsc(dsc_shrink)
save(res,file="res.RData")



