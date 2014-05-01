require("ashr")

#set up some data with mixture of two values of sigmaa
set.seed(100)
sebetahat = 0.01*rgamma(1500,0.5,0.5)
beta =  c(rnorm(500,0,1),rnorm(500,0,0.01),rnorm(500,0,0.000001))
betahat = beta + rnorm(1500,0,sebetahat)

beta.ash = ash(betahat,sebetahat, method="fdr")



zscore = betahat/sebetahat
pval = pchisq(zscore^2,df=1,lower.tail=F)

sum(beta.ash$lfsr<0.05)
#if not installed, first install q value package
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library("qvalue")

qq = qvalue(pval)
sum(qq$qvalues<0.05)



#check whether ordering by q values does better or worse
#job than ordering by confidence, in terms of identifying
#betas with the right sign

err = (sign(betahat) != sign(beta))

plot(cumsum(err[order(qq$qvalues)]),type="l")
lines(cumsum(err[order(beta.ash$lfsr)]),col=2)

#note: I edited nejm_brca_release.txt by removing columns 1-3
hh = read.table("nejm_brca_release_edit.csv",sep=",",skip=3)
subset = apply(hh, 1, max)<20
hh = hh[subset,]

labs = read.table("nejm_brca_release_edit.csv",sep=",",skip=1,nrows=1)
labs = 1*(labs=="BRCA1") + 2 * (labs=="BRCA2") 

hh.betahat = apply(hh[,labs==1],1,mean) - apply(hh[,labs==2],1,mean)
n1 = sum(labs==1)
n2 = sum(labs==2)
hh.sebetahat = sqrt(apply(hh[,labs==1],1,var)/n1 + apply(hh[,labs==2],1,var)/n2)

hh.zscore = hh.betahat/hh.sebetahat
hh.pval = pchisq(hh.zscore^2,df=1,lower.tail=F)
hh.q = qvalue(hh.pval)
sum(hh.q$q<0.05)

