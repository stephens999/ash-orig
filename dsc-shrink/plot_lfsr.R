load("res.RData")
source("set_plot_colors.R")
library(dplyr)



FILE = "../paper/figures/lfsr_est.pdf"

#' @param df dataframe of scores for many methods/scenrios etc
#' @return tall dataframe with columns of scores for each method and the "goldmethod" against which plot is to be made
process_score_for_plotting_against_gold=function(df,PLOTMETHODS=c("ash.n","ash.u","ash.hu"),
                                                 GOLDMETHOD="bayes",PLOTSEEDS=1:100,
                                                 PLOTSCENARIOS=c("spiky","near-normal","flat-top","skew","big-normal","bimodal"),
                                                 PLOTNAMES=PLOTSCENARIOS){
  df %<>% filter(seed %in% PLOTSEEDS) %>% filter(scenario %in% PLOTSCENARIOS) %>% filter(method %in% c(PLOTMETHODS,GOLDMETHOD))
  df$scenario = factor(df$scenario,levels=PLOTSCENARIOS)
  levels(df$scenario)= PLOTNAMES
  #create tall version of dataframe
  df %<>% dplyr::select(-user.self,-sys.self,-elapsed,-user.child,-sys.child) %>%
    reshape2::melt(id.vars=c("method","scenario","seed",".id"),value.name="val")
  #separate bayes and remainder
  df.bayes = df %>% filter(method==GOLDMETHOD)
  df.rest = df %>% filter(method!=GOLDMETHOD)

#join bayes with others, so each line has both the bayes and the non-bayes version
  return(inner_join(df.bayes, df.rest, by=c("scenario","seed","variable")))
}


plot_lfsr=function(lfsr,xlab="True lfsr",ylab="Estimated lfsr",xlim=c(0,0.2),ylim=c(0,0.2),legend.position="bottom"){
  p=ggplot(lfsr,
         aes(val.x,val.y,colour=method.y)) +
    facet_grid(. ~ scenario) + 
    guides(alpha=FALSE) +
    geom_abline(colour = "black") +
    geom_abline(colour= "red", slope=2) +
    xlab(xlab) +
    ylab(ylab) +
    geom_point(shape=1,size=0.1,alpha=0.2) 

  p +scale_y_continuous(limits=ylim) +
        scale_x_continuous(limits=xlim)
         
}

lfsr = process_score_for_plotting_against_gold(res$lfsr,PLOTSEEDS=1:100,PLOTMETHODS="ash.n")
lfdr = process_score_for_plotting_against_gold(res$lfdr,PLOTSEEDS=1:100,PLOTMETHODS="ash.n")

p1=plot_lfsr(lfsr,ylim=c(0,1),xlim=c(0,0.2))
p2=plot_lfsr(lfdr,ylim=c(0,1),xlim=c(0,0.2),xlab="True lfdr",ylab="Estimated lfdr")

# cowplot command: plot_grid(p2+theme_gray()+theme(legend.position="none"),p1+theme_gray()+theme(legend.position="none"),nrow=2,align="v")

print(p1+theme(legend.position="none",axis.text.x = element_text(size = 8,angle=45))
      +coord_equal(ratio=1/5) + colScale)
ggsave("../paper/figures/lfsr_est.png",height=3,width=9)
ggsave("../paper/figures/lfsr_est.pdf",height=3,width=9)


print(p2+theme(legend.position="none",axis.text.x = element_text(size = 8,angle=45))
      +coord_equal(ratio=1/5) + colScale)
ggsave("../paper/figures/lfdr_est.png",height=3,width=9)
ggsave("../paper/figures/lfdr_est.pdf",height=3,width=9)

lfsr.s = process_score_for_plotting_against_gold(res$lfsr,PLOTSEEDS=1:100,PLOTMETHODS="ash.n.s")
p1.s=plot_lfsr(lfsr.s,ylim=c(0,1),xlim=c(0,0.2))

print(p1.s+theme(legend.position="none",axis.text.x = element_text(size = 8,angle=45))
      +coord_equal(ratio=1/5))
ggsave("../paper/figures/lfsr_est_s.png",height=3,width=9)
ggsave("../paper/figures/lfsr_est_s.pdf",height=3,width=9)

lfsr.s.nn = process_score_for_plotting_against_gold(
    res$lfsr,PLOTSEEDS=1:100,PLOTMETHODS="ash.n.s",
    PLOTSCENARIOS=paste0(c("spiky","near-normal","flat-top","skew","big-normal","bimodal"),"-nn"))
p1.s.nn=plot_lfsr(lfsr.s.nn,ylim=c(0,1),xlim=c(0,0.2))

print(p1.s.nn+theme(legend.position="none",axis.text.x = element_text(size = 8,angle=45))
      +coord_equal(ratio=1/5))
ggsave("../paper/figures/lfsr_est_s_nn.png",height=3,width=9)

