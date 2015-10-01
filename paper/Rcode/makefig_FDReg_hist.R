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
qval_fdr = qvalue::lfdr(pval)

pdf("../figures/FDReg_hist.pdf",width=6,height=3)
par(mai=c(0.5,0.5,0.2,0.2),mgp = c(3, 0.5, 0))
layout(matrix(1:2,ncol=2,byrow=TRUE))
plot_FDReg_hist(pval,res.qvalue$pi0,type=1,title="p values",cex.axis=0.8,cex.main=0.8)
plot_FDReg_hist(pval,res.qvalue$pi0,type=4,cex.axis=0.8,yaxt='n',textsize=0.9,cex.main=0.8,title="Decomposition into null/alternative")
axis(side=2, labels=FALSE,tck=-0.01)
dev.off()

############################## decomp_ZA.pdf #############################################

pdf("../figures/decomp_ZA.pdf",width=6,height=6)
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

#a different layout (4 columns, 2 rows) for my Tukey poster
pdf("../figures/decomp_ZA_poster.pdf",width=6,height=2)
layout(matrix(1:8,ncol=4,byrow=TRUE))
par(mar = c(1.3,2,1.5,0.2),mgp = c(3, 0.5, 0))
ncz=25

#p value histograms
altnullhist(pval,qval_fdr,main="p values: qvalue",ncz=ncz,xaxt='n',cex.axis=0.8)
#plot_FDReg_hist(pval,res.qvalue$pi0,type=2,title="p values",xaxt='n',cex.axis=0.8)
axis(side=1, labels=FALSE,tck=-0.01)
#mtext(side=3,"p values",line=1)
altnullhist(pval,res.locfdr$fdr,main="locfdr",ncz=ncz,xaxt='n',cex.axis=0.8)
axis(side=1, labels=FALSE,tck=-0.01)
altnullhist(pval,res.mixfdr$fdr,main="mixfdr",ncz=ncz,xaxt='n',cex.axis=0.8)
axis(side=1, labels=FALSE,tck=-0.01)
altnullhist(pval,res.ash$lfdr,main="ashr",ncz=ncz,xaxt='n',cex.axis=0.8)
axis(side=1, labels=TRUE,tck=-0.01,cex.axis=0.8)
#mtext(side=1,"p values",line=2)

#z score histograms
nullalthist(zscore,qval_fdr,main="z scores: qvalue",ncz=ncz,xaxt='n',cex.axis=0.8,ylim=c(0,0.3),xlim=c(-6,6))
axis(side=1, labels=FALSE,tck=-0.01)
#mtext(side=3,"z scores",line=1)
nullalthist(zscore,res.locfdr$fdr,main="locfdr",ncz=ncz,xaxt='n',cex.axis=0.8,ylim=c(0,0.3),xlim=c(-6,6))
axis(side=1, labels=FALSE,tck=-0.01)
nullalthist(zscore,res.mixfdr$fdr,main="mixfdr",ncz=ncz,xaxt='n',cex.axis=0.8,ylim=c(0,0.3),xlim=c(-6,6))
axis(side=1, labels=FALSE,tck=-0.01)
nullalthist(zscore,res.ash$lfdr,main="ashr",ncz=ncz,xaxt='n',cex.axis=0.8,ylim=c(0,0.3),xlim=c(-6,6))
axis(side=1, labels=TRUE,tck=-0.01,cex.axis=0.8)
#mtext(side=1,"z scores",line=2)

dev.off()




