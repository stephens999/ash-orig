#functions to plot estimates of pi0 from simulations
require(ggplot2)
require(reshape2)

get_pi0.fdrtool = function(f){f$param[3]}
get_pi0.locfdr = function(f){f$fp0[1,3]}
get_pi0.mixfdr = function(f){f$pi[1]}

plot_pi0 = function(sims){
  res=list()
  for(i in 1:length(sims)){
    pi0=sims[[i]]$pi0  
    pi0_ash.n=unlist(lapply(sims[[i]]$betahat.ash.n,get_pi0))
    pi0_ash.u=unlist(lapply(sims[[i]]$betahat.ash.u,get_pi0)) 
    pi0_fdrtool = unlist(lapply(sims[[i]]$betahat.fdrtool,get_pi0.fdrtool))
    pi0_locfdr=unlist(lapply(sims[[i]]$betahat.locfdr,get_pi0.locfdr))
    pi0_mixfdr = unlist(lapply(sims[[i]]$betahat.mixfdr,get_pi0.mixfdr))
    pi0_qval = unlist(lapply(sims[[i]]$betahat.qval,"[[","pi0"))
    
    res[[i]] = data.frame(Scenario=i,pi0=pi0,qvalue=pi0_qval,mixfdr=pi0_mixfdr, locfdr=pi0_locfdr, fdrtool=pi0_fdrtool,ash.nullbiased=pi0_ash.n,ash.uniform=pi0_ash.u)
  }
  cbbPalette <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#000000", "#D55E00", "#CC79A7")
  labels = c('mixfdr','qvalue','locfdr','fdrtool','ash (null-biased)','ash (mle)')
  breaks = c("mixfdr","qvalue","locfdr","fdrtool","ash.nullbiased","ash.uniform")
  
  res.melt = melt(res, id.vars=c("pi0","Scenario"),variable.name="Method")
  res.melt = res.melt[sample(1:nrow(res.melt)),]  
  res.melt$Scenario = as.factor(res.melt$Scenario)
  levels(res.melt$Scenario) = c("Scenario 1a","Scenario 1b","Scenario 2")
  p=ggplot(data=res.melt,aes(pi0,value,colour=Method)) +geom_point(shape=1) +
    facet_grid(. ~ Scenario) +
    geom_abline(colour = "black") +
    xlab("True pi0") +
    ylab("Estimated pi0")
  print(p +scale_y_continuous(limits=c(0,1)) +
          scale_x_continuous(limits=c(0,1)) +
          scale_colour_manual(values=cbbPalette,breaks=breaks,labels=labels) +
          coord_equal(ratio=1))
  
}