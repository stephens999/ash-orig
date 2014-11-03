#This file makes three figures:
#FDReg_hist.pdf: illustrates basic idea of way q value does FDR analysis in p value space
#decomp_ZA.pdf: compares methods in the way they decompose p values and z scores into two groups 
#GOODPOOReg_hist.pdf: illustrates how high-noise measurements can pollute signal


require(ashr)
require(qvalue)
require(locfdr)
require(mixfdr)
source("plot_FDReg_hist.R") #some plotting functions
source("nullalthist.R")

############################ FDReg_hist.pdf ##################################################################

#simple simulated example
ncz = 100 # number of bins in z score histograms
ntest = 10000
set.seed(111)

#simulate (with all tests alternative, true effects $\beta \sim N(0,1)$)
beta = rnorm(ntest)
sebetahat = 1
betahat = beta + rnorm(ntest, 0, sebetahat)
zscore = betahat/sebetahat
pval = pchisq(zscore^2,df=1,lower.tail=F)

#apply the different methods
res.qvalue = qvalue(pval)
res.locfdr = locfdr(zscore,nulltype=0,plot=0)
res.mixfdr = mixFdr(zscore,noiseSD=1,theonull=TRUE,plot=FALSE)
res.ash = ash(betahat,1,method="fdr")

#roughly compute a local fdr for qvalue (to aid the plotting of the decomposition in zscore space)
temp=hist(pval,nclass=50)
bin_fdr = res.qvalue$pi0/temp$density #compute local fdr in each bin of histogram
qval_fdr = bin_fdr[as.numeric(cut(pval,temp$breaks))]
qval_fdr = pmin(1,qval_fdr) #threshold at 1, so fdr never bigger than 1

pdf("figures/FDReg_hist.pdf",width=6,height=3)
par(mai=c(0.5,0.5,0.2,0.2),mgp = c(3, 0.5, 0))
layout(matrix(1:2,ncol=2,byrow=TRUE))
plot_FDReg_hist(pval,res.qvalue$pi0,type=1,title="p values",cex.axis=0.8,cex.main=0.8)
plot_FDReg_hist(pval,res.qvalue$pi0,type=4,cex.axis=0.8,yaxt='n',textsize=0.9,cex.main=0.8,title="Decomposition into null/alternative")
axis(side=2, labels=FALSE,tck=-0.01)
dev.off()

############################## decomp_ZA.pdf #############################################

pdf("figures/decomp_ZA.pdf",width=6,height=6)
layout(matrix(1:12,ncol=3,byrow=FALSE))
plotlabel=function(label,cex=1.5){
  plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  text(0,0,label,cex=cex)
}
plotlabel("qvalue")
plotlabel("locfdr")
plotlabel("mixfdr")
plotlabel("ashr")

par(mar = c(1.3,2,1.5,0.2),mgp = c(3, 0.5, 0))

#p value histograms
altnullhist(pval,qval_fdr,main="p values",ncz=50,xaxt='n',cex.axis=0.8)
#plot_FDReg_hist(pval,res.qvalue$pi0,type=2,title="p values",xaxt='n',cex.axis=0.8)
axis(side=1, labels=FALSE,tck=-0.01)
#mtext(side=3,"p values",line=1)
altnullhist(pval,res.locfdr$fdr,main="",ncz=50,xaxt='n',cex.axis=0.8)
axis(side=1, labels=FALSE,tck=-0.01)
altnullhist(pval,res.mixfdr$fdr,main="",ncz=50,xaxt='n',cex.axis=0.8)
axis(side=1, labels=FALSE,tck=-0.01)
altnullhist(pval,res.ash$lfdr,main="",ncz=50,xaxt='n',cex.axis=0.8)
axis(side=1, labels=TRUE,tck=-0.01,cex.axis=0.8)
#mtext(side=1,"p values",line=2)

#z score histograms
nullalthist(zscore,qval_fdr,main="z scores",ncz=50,xaxt='n',cex.axis=0.8,ylim=c(0,0.3),xlim=c(-6,6))
axis(side=1, labels=FALSE,tck=-0.01)
#mtext(side=3,"z scores",line=1)
nullalthist(zscore,res.locfdr$fdr,main="",ncz=50,xaxt='n',cex.axis=0.8,ylim=c(0,0.3),xlim=c(-6,6))
axis(side=1, labels=FALSE,tck=-0.01)
nullalthist(zscore,res.mixfdr$fdr,main="",ncz=50,xaxt='n',cex.axis=0.8,ylim=c(0,0.3),xlim=c(-6,6))
axis(side=1, labels=FALSE,tck=-0.01)
nullalthist(zscore,res.ash$lfdr,main="",ncz=50,xaxt='n',cex.axis=0.8,ylim=c(0,0.3),xlim=c(-6,6))
axis(side=1, labels=TRUE,tck=-0.01,cex.axis=0.8)
#mtext(side=1,"z scores",line=2)

dev.off()

################################### GOODPOOReg_hist.pdf ###############################

#simple simulated example to illustrate high and low signal
ntest = 10000
set.seed(112)

#simulate (0.5 N(0,1) + 0.5 delta_0)
null_alt = rbinom(ntest,1,0.5)
beta = rnorm(ntest)
beta = ifelse(null_alt==1,beta,0)
GOOD= 1:(ntest/2)
sebetahat = rep(1,ntest)
sebetahat[-GOOD] = 10

betahat = beta + rnorm(ntest, 0, sebetahat)
zscore = betahat/sebetahat
pval = pchisq(zscore^2,df=1,lower.tail=F)

pdf("figures/GOODPOOReg_hist.pdf",width=6.5,height=3)
par(mai=c(0.3,0.3,0.2,0.2),mgp = c(3, 0.5, 0))
layout(matrix(1:3,ncol=3,byrow=TRUE))
plot_FDReg_hist(pval[GOOD],1,type=1,title="Low-noise observations",ylab="",nc=20,cex.axis=1,cex.main=1.2,ylim=c(0,2.5))
plot_FDReg_hist(pval[-GOOD],1,type=1,title="High-noise observations",ylab="",nc=20,yaxt='n',cex.axis=1,cex.main=1.2,ylim=c(0,2.5))
axis(side=2, labels=FALSE,tck=-0.01)
plot_FDReg_hist(pval,1,type=1,title="Combined",yaxt='n',ylab="",nc=20,cex.axis=1,cex.main=1.2,ylim=c(0,2.5))
axis(side=2, labels=FALSE,tck=-0.01)
dev.off()

################################### GOODPOOReg_scatterplot.pdf ###############################

#apply the different methods
res.qvalue = qvalue(pval)
res.locfdr = locfdr(zscore,nulltype=0,plot=0)
res.ash = ash(betahat,sebetahat,method="fdr")

res.qvalue.good = qvalue(pval[GOOD])
res.locfdr.good = locfdr(zscore[GOOD],nulltype=0,plot=0)
res.ash.good = ash(betahat[GOOD],sebetahat[GOOD],method="fdr")

pdf("figures/GOODPOOReg_scatter.pdf",width=6.5,height=3)
par(mai=c(0.3,0.3,0.2,0.2),mgp = c(3, 0.5, 0))
layout(matrix(1:3,ncol=3,byrow=TRUE))
plot(res.qvalue.good$q,res.qvalue$q[GOOD],main="qvalue",xlim=c(0,1),ylim=c(0,1),axes=F)
axis(side=2)
axis(side=1)
abline(a=0,b=1,col=2)
plot(res.locfdr.good$fdr,res.locfdr$fdr[GOOD],main="locfdr",xlim=c(0,1),ylim=c(0,1),axes=F)
axis(side=1)
abline(a=0,b=1,col=2)
plot(res.ash.good$lfsr,res.ash$lfsr[GOOD],main="ash",xlim=c(0,1),ylim=c(0,1),axes=F)
axis(side=1)
abline(a=0,b=1,col=2)
dev.off()

res = rbind(data.frame(x=res.ash.good$lfsr,y=res.ash$lfsr[GOOD],type="ash"), 
            data.frame(x=res.locfdr.good$fdr,y=res.locfdr$fdr[GOOD],type='locfdr'), 
            data.frame(x= res.qvalue.good$qvalues, y = res.qvalue$qvalues[GOOD],type="qvalue") )

library("ggplot2")
pdf("figures/GOODPOOReg_scatter.pdf",height=3,width=6.5)
pp= ggplot(data=res,aes(x,y)) +geom_point(shape=1) +
  facet_grid(. ~ type) +
  geom_abline(colour = "red") +
  xlab("Analysing low-noise data only") +
  ylab("Analysing combined data")


print(pp +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1))  +
        coord_equal(ratio=1))

dev.off()
