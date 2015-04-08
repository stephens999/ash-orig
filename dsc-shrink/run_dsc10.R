library(dscr)

dsc_shrink=new.dsc("shrink10","dsc-shrink10-files")
source("addScenarios10.R")
source("addMethods.R")

addOutputParser(dsc_shrink,"ash2beta",ash2beta_est,"ash_output","beta_est_output")
addOutputParser(dsc_shrink,"mixfdr2beta",mixfdr2beta_est,"mixfdr_output","beta_est_output")

addOutputParser(dsc_shrink,"ash2pi0",ash2pi0_est,"ash_output","pi0_est_output")
addOutputParser(dsc_shrink,"mixfdr2pi0",mixfdr2pi0_est,"mixfdr_output","pi0_est_output")
addOutputParser(dsc_shrink,"locfdr2pi0",locfdr2pi0_est,"locfdr_output","pi0_est_output")
addOutputParser(dsc_shrink,"qvalue2pi0",qvalue2pi0_est,"qvalue_output","pi0_est_output")

addOutputParser(dsc_shrink,"ash2fitted.g",ash2fitted.g,"ash_output","g_output")
#addOutputParser(dsc_shrink,"mixfdr2fitted.g",mixfdr2fitted.g,"mixfdr_output","g_output")



source("score.R")
addScore(dsc_shrink,score,"beta_err","beta_est_output")
addScore(dsc_shrink,score2,"pi0_score","pi0_est_output")
addScore(dsc_shrink,score3,"cdf_score","g_output")
addScore(dsc_shrink,score_neg,"negprob","ash_output") #just extract the negativeprobs
addScore(dsc_shrink,score_pos,"posprob","ash_output") #just extracts the positiveprobs
addScore(dsc_shrink,score_fdr,"fdr","mixfdr_output") #just extracts the fdr
addScore(dsc_shrink,score_betahat,"betahat","mixfdr_output") #just extracts the fdr



res10=run_dsc(dsc_shrink)
save(res10,dsc_shrink,file="res10.RData")



