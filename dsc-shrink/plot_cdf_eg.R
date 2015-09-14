library(plyr)
library(dplyr)
library(dscr)
library(magrittr)
library(ashr)
library(ggplot2)
source("methods/mixfdr.wrapper.R")
source("set_plot_colors.R")
load("res.RData")

#df is a list with components method,seed, scenario
#cdf is evaluated at x
#returns list of x,y,pi0, where cdf values are in y
get_cdf = function(df,dsc=dsc_shrink,x=seq(-6,6,length=100)){
  m=df$method
  if(m=="truth"){m="ash.n"}
  temp=loadExample(dsc,df$seed,df$scenario,m)
  pi0=temp$meta$pi0
  
  if(df$method=="truth"){
    s=dsc$scenarios[[df$scenario]]
    galt=s$args$g
    g=normalmix(c(pi0,(1-pi0)*mixprop(galt)),
                c(0,comp_mean(galt)),
                c(0,comp_sd(galt)))
  } else {
    if(grepl("mixfdr",df$method)){
      temp$output= mixfdr2fitted.g(temp$output)
    }
    g=temp$output$fitted.g
  }
  data.frame(x=x,y=as.numeric(mixcdf(g,x)),pi0=pi0)
}

plot_mean_cdf = function(SEEDS,PLOTMETHODS=c("ash.n","ash.u","ash.hu","truth","mixfdr.tnull"),
      PLOTSCENARIOS=c("spiky","near-normal","flat-top","skew","big-normal","bimodal"),pi0filter=FALSE,...){
  PLOTNAMES=PLOTSCENARIOS
  #set up dataframe with cdf for all methods and all datasets
  df= expand.grid(seed=SEEDS,scenario=PLOTSCENARIOS, method=PLOTMETHODS,stringsAsFactors = FALSE)
  df.cdf=ddply(df,.(seed,scenario,method), get_cdf)
  
  if(pi0filter==TRUE){
    df.cdf %<>% filter(pi0<0.55 & pi0>0.45)
  }
  if(length(SEEDS)>1){
    df.cdf %<>% group_by(x,method,scenario) %>% dplyr::summarise(y=mean(y))
  }
  
  df.cdf$scenario=factor(df.cdf$scenario,levels=PLOTSCENARIOS)
  levels(df.cdf$scenario)=PLOTNAMES

  ggplot(df.cdf,aes(x=x,y=y,color=method),...) + colScale + geom_line() + facet_grid(.~scenario) + theme(legend.position = "bottom")

}

pdf("../paper/figures/egcdf.pdf",height=3,width=9)
plot_mean_cdf(1)
dev.off()

pdf("../paper/figures/mean_cdf.pdf",height=3,width=9)
plot_mean_cdf(1:100,PLOTMETHODS=c("ash.n","ash.u","ash.hu","truth"),pi0filter=TRUE)
dev.off()

pdf("../paper/figures/mean_cdf_nopen.pdf",height=3,width=9)
names(myColors) <- c("mixfdr.tnull","ash.hu.s","ash.n.s","ash.u.s","qvalue","locfdr","truth")
colScale <- scale_colour_manual(name = "method",values = myColors)
plot_mean_cdf(1:100,PLOTMETHODS=c("ash.n.s","ash.u.s","ash.hu.s","truth"),pi0filter=TRUE)
dev.off()
