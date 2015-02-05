require(ashr)
require(qvalue)
require(fdrtool)
require(mixfdr)
require(locfdr)
require(ggplot2)   

load("sim1.RData")
load("simABC.RData")
source("plot_examples.R")
source("plot_rmse.R")
source("plot_pi0.R")
source("plot_lfsr.R")
source("plot_FDReg_hist.R")



pdf("figures/simABC_egdens.pdf",width=6.5,height=2)
plot_examples(list(simA,simB,simC))
dev.off()

pdf("figures/rmse_biplot.pdf",width=6.5,height=2)
plot_rmse(list(simA,simB,simC))
dev.off()

#This figure not used in paper?
pdf("figures/rmse_biplot_withzerobetahat.pdf")
plot_rmse(list(simA,simB,simC),inczero=TRUE,incbetahat=TRUE)
dev.off()

pdf("figures/rmse_boxplot.pdf",width=6.5,height=2)
plot_rmse_boxplot(list(simA,simB,simC))
dev.off()

pdf("figures/rmse_boxplot_extended0.pdf",width=6.5,height=2)
plot_rmse_boxplot(list(simA,simB,simC),TRUE,TRUE,TRUE)
dev.off()

pdf("figures/rmse_boxplot_extended1.pdf",width=6.5,height=2)
plot_rmse_boxplot(list(simA,simB,simC),incmixfdr=TRUE)
dev.off()


pdf("figures/loglik_boxplot.pdf",width=6.5,height=2)
plot_loglik_boxplot(list(simA,simB,simC))
dev.off()

pdf("figures/rmse_loglik_boxplot.pdf",width=6.5,height=5)
plot_rmse_loglik_boxplot(list(simA,simB,simC))
dev.off()


pdf("figures/simABC_eg_withfit.pdf",width=6.5,height=2)
plot_examples_withfit(list(simA,simB,simC))
dev.off()

pdf("figures/simABC_eg_cdf_withfit.pdf",width=6.5,height=2)
plot_examples_cdf_withfit(list(simA,simB,simC))
dev.off()

plot_LR=function(sims){
  hist(unlist(lapply(sims$fit.ash.u,get_loglik))-unlist(lapply(sims$fit.ash.n,get_loglik)), xlab="loglik difference", main="loglik differences for nullbiased prior vs mle",nclass=10)
}

pdf("figures/logLR.pdf")
plot_LR(simA)
plot_LR(simB)
dev.off()




#Plot pi0 from each method


png("figures/estpi0_sim1sim2.png",height=160,width=540)
plot_pi0(list(simres1a,simres1b,simres2))
dev.off()



png("figures/lfsdr_sim1sim2_blowup.png",height=427,width=720)
plot_lfsdr(list(simres1a,simres1b,simres2),0.1,ptype="lfsra")
dev.off()

png("figures/lfsdr_sim1sim2_blowup.png",height=427,width=720)
plot_lfsdr(list(simres1a,simres1b,simres2),0.1,ymax=0.2,ptype="lfsra")
dev.off()

png("figures/lfdr_sim1sim2_blowup.png",height=160,width=540)
plot_lfsr(list(simres1a,simres1b,simres2),0.1,ptype="lfdr")
dev.off()

png("figures/lfsra_sim1sim2_blowup.png",height=160,width=540)
plot_lfsr(list(simres1a,simres1b,simres2),0.1,ptype="lfsra")
dev.off()

png("figures/lfsra_sim1sim2_blowup.png",height=160,width=540)
plot_lfsr(list(simres1a,simres1b,simres2),0.1,ptype="lfsr")
dev.off()


