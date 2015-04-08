library(plyr)
library(dplyr)
library(dscr)
library(magrittr)
library(ashr)
library(ggplot2)
source("methods/mixfdr.wrapper.R")
load("res.RData")

MEANPLOTFILE = "../paper/figures/meancdf.pdf"
PLOTSCENARIOS = c("hard","An","Bn","Cn","easy")
PLOTNAMES = c("spiky","near-normal","flat-top","skew","big-normal")
PLOTMETHODS = c("ash.n","ash.hu","truth","mixfdr.tnull")
SEEDS = 1

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

#set up dataframe with cdf for all methods and all datasets
df= expand.grid(seed=SEEDS,scenario=PLOTSCENARIOS, method=PLOTMETHODS,stringsAsFactors = FALSE)

df.cdf=ddply(df,.(seed,scenario,method), get_cdf)
df.cdf.mean = df.cdf %>% group_by(x,method,scenario) %>% summarise(y=mean(y))
df.midpi0.cdf.mean = df.cdf %>% filter(pi0<0.55 & pi0>0.45) %>% group_by(x,method,scenario) %>% summarise(y=mean(y))

pdf(MEANPLOTFILE)

df.cdf$scenario=factor(df.cdf$scenario,levels=PLOTSCENARIOS)
levels(df.cdf$scenario)=PLOTNAMES

df.cdf$seed = factor(df.cdf$seed)

ggplot(df.cdf,
       aes(x=x,y=y,color=method)) + geom_line(aes(linetype=seed)) + facet_grid(.~scenario)

ggplot(df.cdf.mean,
       aes(x=x,y=y,color=method)) + geom_line() + facet_grid(.~scenario)


ggplot(df.midpi0.cdf.mean,
       aes(x=x,y=y,color=method)) + geom_line() + facet_grid(.~scenario)



dev.off()

