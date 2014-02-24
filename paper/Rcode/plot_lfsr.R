require(reshape2)
require(ggplot2)

mean_quant=function(x, mult = 1){
  x <- na.omit(x)
  sd <- mult * sqrt(var(x))
  mean <- mean(x)
  data.frame(y = median(x), ymin = quantile(x,0.25) , ymax = quantile(x,0.75))
}

#ptype indicates what type of plot to do
#maxlfsr controls maximum x axis value
#maxy controls maximum y axis value
plot_lfsr=function(sims,maxlfsr=0.1,ptype=c("lfsr","lfsra","lfdr"),maxy=1){
  ptype=match.arg(ptype)
  xlabtype = ifelse(ptype=="lfdr","lfdr","lfsr")
  res=list()
  for(i in 1:length(sims)){
    lfsr.ash.n= unlist(lapply(sims[[i]]$betahat.ash.n,"[[",ptype))
    lfsr.ash.u= unlist(lapply(sims[[i]]$betahat.ash.u,"[[",ptype))
    if(ptype=="lfdr"){
      lfsr.ash.true= unlist(lapply(sims[[i]]$betahat.ash.true,"[[","lfdr"))
      lfdr.mixfdr = unlist(lapply(sims[[i]]$betahat.mixfdr,"[[","fdr"))
      lfdr.locfdr = unlist(lapply(sims[[i]]$betahat.locfdr,"[[","fdr"))
    } else {
      lfsr.ash.true= unlist(lapply(sims[[i]]$betahat.ash.true,"[[","lfsr"))
    }
    
    subset = lfsr.ash.true<maxlfsr
    
    if(length(subset)>0){  
      res[[i]] = data.frame(Scenario=i,ash.nullbiased=lfsr.ash.n[subset],
                            ash.uniform=lfsr.ash.u[subset],
                            Bayes=0.1*maxlfsr*findInterval(lfsr.ash.true[subset],seq(0,maxlfsr,length=11))-0.05*maxlfsr)
      if(ptype=="lfdr"){
        res[[i]]=data.frame(res[[i]],mixfdr = lfdr.mixfdr[subset])
      }
    }
  }
  

  
  res.melt = melt(res, id.vars=c("Bayes","Scenario"),variable.name="Method")
  
  
  cbbPalette <- c("#000000","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  labels = c('ash (null-biased)','ash (mle)')
  breaks = c("ash.nullbiased","ash.uniform")
  if(ptype=="lfdr"){
    labels = c('mixfdr', labels)
    breaks = c("mixfdr", breaks)
  }
  
  p=ggplot(data=res.melt,aes(Bayes,value,colour=Method)) + 
    facet_grid(. ~ Scenario) +
    #scale_fill_manual(values=cbbPalette) +
    #scale_colour_manual(values=cbbPalette) +
    geom_point(size=1,alpha=0.1) + 
    stat_smooth(se=FALSE,size=2) +
    stat_summary(fun.data="mean_quant", geom="ribbon", alpha=0.25) +
    geom_abline(colour = "black") +
    xlab(paste0(xlabtype," (Bayes)")) +
    ylab(paste0(ptype," (Method)")) 
  
  print(p +scale_y_continuous(limits=c(0,maxy)) +
          scale_x_continuous(limits=c(0,maxlfsr)) +
          scale_colour_manual(values=cbbPalette,breaks=breaks,labels=labels))         
}


#plots estimated lfdr and lfsr against truth.
#ptype indicates whether to use lfsr or lfsra
#maxx controls maximum x axis value
#maxy controls maximum y axis value
plot_lfsdr=function(sims,maxx=0.1,ptype=c("lfsr","lfsra"),maxy=1){
  ptype=match.arg(ptype)
  res=list()
  res.lfsr=list()
  res.lfdr=list()
  
  for(i in 1:length(sims)){
    lfsr.ash.true= unlist(lapply(sims[[i]]$betahat.ash.true,"[[","lfsr"))
    lfsr.ash.n= unlist(lapply(sims[[i]]$betahat.ash.n,"[[",ptype))
    lfsr.ash.u= unlist(lapply(sims[[i]]$betahat.ash.u,"[[",ptype))
    
    lfdr.ash.n= unlist(lapply(sims[[i]]$betahat.ash.n,"[[","lfdr"))
    lfdr.ash.u= unlist(lapply(sims[[i]]$betahat.ash.u,"[[","lfdr"))
    lfdr.ash.true= unlist(lapply(sims[[i]]$betahat.ash.true,"[[","lfdr"))
    lfdr.mixfdr = unlist(lapply(sims[[i]]$betahat.mixfdr,"[[","fdr"))
    lfdr.locfdr = unlist(lapply(sims[[i]]$betahat.locfdr,"[[","fdr"))
    
    
    subset = lfsr.ash.true<maxx
    
    res.lfsr[[i]] = data.frame(Scenario=i,Measure='lfsr', ash.nullbiased=lfsr.ash.n[subset], ash.uniform=lfsr.ash.u[subset],Bayes=0.1*maxx*findInterval(lfsr.ash.true[subset],seq(0,maxx,length=11))-0.05*maxx,mixfdr = NA)
    
    subset = lfdr.ash.true<maxx
    res.lfdr[[i]] = data.frame(Scenario=i,Measure='lfdr', ash.nullbiased=lfdr.ash.n[subset], ash.uniform=lfdr.ash.u[subset],Bayes=0.1*maxx*findInterval(lfdr.ash.true[subset],seq(0,maxx,length=11))-0.05*maxx,mixfdr = lfdr.mixfdr[subset])
    
    
    res[[i]]=rbind(res.lfdr[[i]],res.lfsr[[i]])
  }
  
  res.melt = melt(res, id.vars=c("Bayes","Scenario","Measure"),variable.name="Method")
  res.melt$Scenario = as.factor(res.melt$Scenario)
  levels(res.melt$Scenario) = c("Scenario 1a","Scenario 1b","Scenario 2")
  
  cbbPalette <- c("#E69F00","#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  labels = c('mixfdr','ash (null-biased)','ash (mle)')
  breaks = c("mixfdr","ash.nullbiased","ash.uniform")
  
  p=ggplot(data=res.melt,aes(Bayes,value,colour=Method)) + 
    facet_grid(Measure ~ Scenario ) +
    #scale_fill_manual(values=cbbPalette) +
    #scale_colour_manual(values=cbbPalette) +
    geom_point(size=1,alpha=0.1) + 
    #stat_smooth(se=FALSE,size=2) +
    stat_summary(fun.data="mean_quant", geom="ribbon", alpha=0.25) +
    geom_abline(colour = "red",size=1) +
    xlab("Truth") +
    ylab("Estimate") 
  
  print(p +scale_y_continuous(limits=c(0,maxy)) +
          scale_x_continuous(limits=c(0,maxx)) +
          scale_colour_manual(values=cbbPalette,breaks=breaks,labels=labels))        
}
