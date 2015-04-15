load("res.RData")
source("set_plot_colors.R")
library(dplyr)

PLOTMETHODS= c("mixfdr.tnull","ash.n","ash.u","qvalue","locfdr")
PLOTSCENARIOS = c("spiky","near-normal","flat-top","skew","big-normal","bimodal")
PLOTNAMES = PLOTSCENARIOS
ALPHALEVEL = 0.8 # controls transparency
FILE = "../paper/figures/pi0_est.pdf"


pdf(FILE,height=3,width=9)

df = res$pi0_score %>% filter(scenario %in% PLOTSCENARIOS) %>% filter(method %in% PLOTMETHODS)
df$scenario = factor(df$scenario,levels=PLOTSCENARIOS)
levels(df$scenario)= PLOTNAMES

p=ggplot(df,
         aes(pi0,pi0_est,colour=method,alpha=ALPHALEVEL)) +geom_point(shape=1) +
  facet_grid(. ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(colour = "black") +
  xlab("True pi0") +
  ylab("Estimated pi0") 
print(p +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        #scale_colour_manual(values=cbbPalette,breaks=breaks,labels=labels) +
        coord_equal(ratio=1) + colScale)

dev.off()
