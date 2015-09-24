#the aim of this dsc is to assess the distribution of the logLR 
#value for ash under the null that all effects are 0

library(dscr)

dsc_znull=new.dsc("znull","dsc-znull-files")
source("add_scenarios.null.R")
source("add_methods.null.R")

source("score.R")
add_score(dsc_znull,score_logLR,"logLR","ash_output")


res=run_dsc(dsc_znull)
save(res,dsc_znull,file="res.znull.RData")



library(ggplot2)
library(dplyr)
temp=ungroup(res %>% group_by(scenario,method) %>% summarise(gt0=mean(logLR>0)))
xtabs(gt0~method+scenario,temp)

# scenario
# method            znull.100 znull.1000
# ash.hu              0.438      0.742
# ash.hu.nw2          0.814      0.916
# ash.hu.s            0.988      0.998
# ash.hu.s.gmfine     0.978      0.996
# ash.n               0.170      0.352
# ash.n.nw2           0.330      0.496
# ash.n.s             0.476      0.490
# ash.n.s.gmfine      0.472      0.502
# ash.u               0.210      0.424
# ash.u.nw2           0.356      0.534
# ash.u.s             0.482      0.570
# ash.u.s.gmfine      0.484      0.562
# 
# Note that the ash.n.s basically conforms exactly to the null expectation
# of 0.5 chi2_0 +0.5 chi2_1 from Stram and Lee, Biometrics
# conjecture: the uniform has the same asymptotic behaviour
# and the half uniform is like the sum of two of these?

#from this
#> qchisq(0.9,df=1)
#[1] 2.705543
# so a 95% procedure would need to check for logLR > 2.7/2 = 1.35


pdf("qplots.znull.pdf")
ggplot(res %>% filter(logLR!=0) %>% filter(grepl("ash.hu",method)) , aes(sample = logLR)) + facet_grid(scenario~method) +
  stat_qq(dist = qchisq, dparams=list(df=1)) + geom_abline() + ggtitle("qqplot under null vs chisq-1; line slope=1")
ggplot(res %>% filter(logLR!=0) %>% filter(grepl("ash.u",method)) , aes(sample = logLR)) + facet_grid(scenario~method) +
  stat_qq(dist = qchisq, dparams=list(df=1)) + geom_abline(intercept=0,slope=0.5) + ggtitle("qqplot under null vs chisq-1; line slope=0.5")
ggplot(res %>% filter(logLR!=0) %>% filter(grepl("ash.n",method)), aes(sample = logLR)) + facet_grid(scenario~method) +
  stat_qq(dist = qchisq, dparams=list(df=1)) + geom_abline(intercept=0,slope=0.5) + ggtitle("qqplot under null vs chisq-1; line slope=0.5")
dev.off()



