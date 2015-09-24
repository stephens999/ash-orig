library(dscr)

dsc_shrink=new.dsc("shrink100","dsc-shrink100-files")
source("add_scenarios100.R")
source("add_methods.R")

add_output_parser(dsc_shrink,"ash2beta",ash2beta_est,"ash_output","beta_est_output")
add_output_parser(dsc_shrink,"mixfdr2beta",mixfdr2beta_est,"mixfdr_output","beta_est_output")

add_output_parser(dsc_shrink,"ash2pi0",ash2pi0_est,"ash_output","pi0_est_output")
add_output_parser(dsc_shrink,"mixfdr2pi0",mixfdr2pi0_est,"mixfdr_output","pi0_est_output")

add_output_parser(dsc_shrink,"ash2fitted.g",ash2fitted.g,"ash_output","g_output")
#add_output_parser(dsc_shrink,"mixfdr2fitted.g",mixfdr2fitted.g,"mixfdr_output","g_output")



source("score.R")
add_score(dsc_shrink,score,"beta_err","beta_est_output")
add_score(dsc_shrink,score2,"pi0_score","pi0_est_output")
add_score(dsc_shrink,score3,"cdf_score","g_output")
add_score(dsc_shrink,score_neg,"negprob","ash_output") #just extract the negativeprobs
add_score(dsc_shrink,score_pos,"posprob","ash_output") #just extracts the positiveprobs
add_score(dsc_shrink,score_fdr,"fdr","mixfdr_output") #just extracts the fdr
add_score(dsc_shrink,score_betahat,"betahat","mixfdr_output") #just extracts the fdr



res100=run_dsc(dsc_shrink)
save(res100,dsc_shrink,file="res100.RData")


xtabs(pi0_est ~ method+scenario,res100$pi0_score)
xtabs(user.self ~ method+scenario,res100$pi0_score)


