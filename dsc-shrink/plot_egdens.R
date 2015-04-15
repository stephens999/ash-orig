library("ashr")
library("ggplot2")
library("dplyr")
library("dscr")
load("res.RData")
PLOTSCENARIOS = c("spiky","near-normal","flat-top","skew","big-normal","bimodal")
PLOTNAMES = PLOTSCENARIOS
PLOTFILE = "../paper/figures/scenarios_density.pdf"

df=data.frame()

for(i in PLOTSCENARIOS){
  s=dsc_shrink$scenarios[[i]]
  g=s$args$g
  x = seq(-6,6,length=100)
  y = as.numeric(dens(g,x))
  df = rbind(df,data.frame(x=x,y=y,scenario=i))
}


df$scenario = factor(df$scenario,levels=PLOTSCENARIOS)
levels(df$scenario)= PLOTNAMES
  
pdf(PLOTFILE,height=3,width=9)
ggplot(df, aes(x=x,y=y)) + geom_line(size=1.2,linetype=1) + facet_grid(.~scenario) + ylab("density")
dev.off()


  
 
