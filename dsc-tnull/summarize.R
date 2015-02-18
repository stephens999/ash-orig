library(ggplot2)
library(dplyr)

res$err = res$pi0_est<0.95
res.sum1=aggregate(err~method+scenario,res, mean)

xtabs(err ~ method+scenario,res.sum1)
xtabs(fp1 ~ method+scenario,res)
xtabs(fp2 ~ method+scenario,res)

