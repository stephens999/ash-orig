require(ggplot2)

cdf.mixfdr = function(a,x){mixcdf(normalmix(a$pi,a$mu,a$sigma-1),x)}

plot_examples = function(sims){
  sim1= sims[[1]]
  sim2= sims[[2]]
  sim3= sims[[3]]
  len=length(sim1$beta[[1]])
  res = data.frame(beta = c(sim1$beta[[1]],sim2$beta[[1]],sim3$beta[[1]]), Scenario = c(rep(1,len),rep(2,len),rep(3,len)))
  x = seq(-8,6,length=100)
  y = dnorm(x)
  dens1 = data.frame(x=x,dtrue=t(density(sim1$fit.ash.true[[1]],x)$y), Scenario=1)
  dens2 = data.frame(x=x,dtrue=t(density(sim2$fit.ash.true[[1]],x)$y), Scenario=2)
  dens3 = data.frame(x=x,dtrue=t(density(sim3$fit.ash.true[[1]],x)$y), Scenario=3)
  dens=rbind(dens1,dens2,dens3)
  ggplot(res) + facet_grid(.~Scenario)  + 
    geom_line(data = dens, aes(x,dtrue),size=1.2,alpha=0.9,linetype=1) + scale_x_continuous(limits=c(-6,6))
}

plot_examples_cdf_withfit = function(sims){
  sim1= sims[[1]]
  sim2= sims[[2]]
  sim3= sims[[3]]
  len=length(sim1$beta[[1]])
  res = data.frame(beta = c(sim1$beta[[1]],sim2$beta[[1]],sim3$beta[[1]]), Scenario = c(rep(1,len),rep(2,len),rep(3,len)))
  x = seq(-6,6,length=100)
  y = dnorm(x)
  cdf1 = data.frame(x=x,dn=t(cdf.ash(sim1$fit.ash.n[[1]],x)$y),du=t(cdf.ash(sim1$fit.ash.u[[1]],x)$y),dhu=t(cdf.ash(sim1$fit.ash.hu[[1]],x)$y), dtrue=t(cdf.ash(sim1$fit.ash.true[[1]],x)$y), dn.fdr=t(cdf.ash(sim1$fit.ash.fdr.n[[1]],x)$y),dm = t(cdf.mixfdr(sim1$fit.mixfdr[[1]],x)), Scenario=1)
  cdf2 = data.frame(x=x,dn=t(cdf.ash(sim2$fit.ash.n[[1]],x)$y),du=t(cdf.ash(sim2$fit.ash.u[[1]],x)$y),dhu=t(cdf.ash(sim2$fit.ash.hu[[1]],x)$y), dtrue=t(cdf.ash(sim2$fit.ash.true[[1]],x)$y),dn.fdr=t(cdf.ash(sim2$fit.ash.fdr.n[[1]],x)$y),dm = t(cdf.mixfdr(sim2$fit.mixfdr[[1]],x)),Scenario=2)
  cdf3 = data.frame(x=x,dn=t(cdf.ash(sim3$fit.ash.n[[1]],x)$y),du=t(cdf.ash(sim3$fit.ash.u[[1]],x)$y),dhu=t(cdf.ash(sim3$fit.ash.hu[[1]],x)$y), dtrue=t(cdf.ash(sim3$fit.ash.true[[1]],x)$y),dn.fdr=t(cdf.ash(sim3$fit.ash.fdr.n[[1]],x)$y),dm = t(cdf.mixfdr(sim3$fit.mixfdr[[1]],x)),Scenario=3)
  
  cdf=rbind(cdf1,cdf2,cdf3)
  cbbPalette <- c("#56B4E9", "#009E73","#F0E442",  "#000000","#0072B2", "#E69F00", "#000000", "#D55E00", "#CC79A7")
  ggplot(res, aes(x=beta)) + facet_grid(.~Scenario)  +
    geom_line(data = cdf, aes(x,dtrue,color='True'),size=1.5,alpha=0.9,linetype=1) +
    geom_line(data = cdf, aes(x,du,color='Uniform'),size=1,alpha=0.9,linetype=2) +
    geom_line(data = cdf, aes(x,dhu,color='Half Uniform'),size=1,alpha=0.9,linetype=2) +
    geom_line(data = cdf, aes(x,dn,color='Normal'),size=1,alpha=0.9,linetype=2) +
    geom_line(data = cdf, aes(x,dn.fdr,color='Normal, null-biased'),size=1,alpha=0.9,linetype=2) +
    geom_line(data = cdf, aes(x,dm,color='mixfdr'),size=1,alpha=0.9,linetype=2) +
    scale_colour_manual(name = 'Method', values = cbbPalette) + scale_x_continuous(limits=c(-6,6))
}


plot_examples_withfit = function(sims){
  sim1= sims[[1]]
  sim2= sims[[2]]
  sim3= sims[[3]]
  len=length(sim1$beta[[1]])
  res = data.frame(beta = c(sim1$beta[[1]],sim2$beta[[1]],sim3$beta[[1]]), Scenario = c(rep(1,len),rep(2,len),rep(3,len)))
  x = seq(-8,6,length=100)
  y = dnorm(x)
  dens1 = data.frame(x=x,dn= t(density(sim1$fit.ash.n[[1]],x)$y),du=t(density(sim1$fit.ash.u[[1]],x)$y),dhu=t(density(sim1$fit.ash.hu[[1]],x)$y), Scenario=1)
  dens2 = data.frame(x=x,dn = t(density(sim2$fit.ash.n[[1]],x)$y),du=t(density(sim2$fit.ash.u[[1]],x)$y), dhu=t(density(sim2$fit.ash.hu[[1]],x)$y), Scenario=2)
  dens3 = data.frame(x=x,dn = t(density(sim3$fit.ash.n[[1]],x)$y),du=t(density(sim3$fit.ash.u[[1]],x)$y), dhu=t(density(sim3$fit.ash.hu[[1]],x)$y), Scenario=3)
  dens=rbind(dens1,dens2,dens3)
  cbbPalette <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#000000", "#D55E00", "#CC79A7")
  ggplot(res, aes(x=beta)) + facet_grid(.~Scenario)  +
    geom_histogram(aes(y=..density..),binwidth=0.5,
                   colour="black", fill="white") +  
    geom_line(data = dens, aes(x,du,color='Uniform'),size=1.2,alpha=0.9,linetype=1) +
    geom_line(data = dens, aes(x,dhu,color='Half Uniform'),size=1.2,alpha=0.9,linetype=1) +
    geom_line(data = dens, aes(x,dn,color='Normal'),size=1.2,alpha=0.9,linetype=1) +
    scale_colour_manual(name = 'Method', values = cbbPalette) + scale_x_continuous(limits=c(-4,4))
}


