library("ashr")
require(reshape2)
require(ggplot2)


##Global parameters for generating the data to feed ash() function
mixsds = c(0,0.25,0.5,1,2)
mixpi_alt = c(0.4,0.2,0.2,0.2) #mixture proportions under the alternative
set.seed(2014)


##functions for simulation and generating plots
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function()
{
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
}

unimix = function(pi,a,b){
  structure(data.frame(pi,a,b),class="unimix")
}

normalmix = function(pi,mean,sd){
  structure(data.frame(pi,mean,sd),class="normalmix")
}

nonzerosim=function(mixsd,mixpi_alt,betamean=5,bsd=1,method="shrink",mixcompdist="normal",minpi0=0,seedval=2014,nsamp=1000,niter=20){
  set.seed(seedval)
  beta =list()
  betahatsd=list()
  betahat = list()
  
  betahat.ash.nonzero = list()
  betahat.ash.zero = list()
  betahat.ash.true=list() 
  pi0 = rep(0,niter)	
  
  for(i in 1:niter){
    pi0[i]=runif(1,minpi0,1)
    mixpi = c(pi0[i],(1-pi0[i])*mixpi_alt)
    sds = sample(mixsd,nsamp,prob=mixpi,replace=TRUE)
    beta[[i]] = rnorm(nsamp,betamean,sds)
    betahatsd[[i]] = bsd
    betahat[[i]] = beta[[i]]+rnorm(nsamp,0,betahatsd[[i]])
    
    betahat.ash.zero[[i]] = ash(betahat[[i]],betahatsd[[i]],method=method,mixcompdist=mixcompdist,nonzeromean=FALSE)
    betahat.ash.nonzero[[i]] = ash(betahat[[i]],betahatsd[[i]],method=method,mixcompdist=mixcompdist,nonzeromean=TRUE)
    betahat.ash.true[[i]] = ash(betahat[[i]]-betamean,betahatsd[[i]],nonzeromean=FALSE,g=normalmix(mixpi,rep(0,length(mixpi)),mixsd))
    betahat.ash.true[[i]]$PosteriorMean=betahat.ash.true[[i]]$PosteriorMean+betamean
	if(class(betahat.ash.true[[i]]$fitted.g)=="normalmix"){
		betahat.ash.true[[i]]$fitted.g$mean=rep(betamean,length(betahat.ash.true[[i]]$fitted.g$pi))
	}else{
		betahat.ash.true[[i]]$fitted.g$a=betahat.ash.true[[i]]$fitted.g$a+betamean
		betahat.ash.true[[i]]$fitted.g$b=betahat.ash.true[[i]]$fitted.g$b+betamean
	}
    print(betahat.ash.true[[i]]$PosteriorMean[1:20])
    print(i)
    toc()
}
  return(list(beta =beta,
              betahatsd=betahatsd,
              betahat = betahat,        
			  betahat.ash.nonzero=  betahat.ash.nonzero,
     		  betahat.ash.zero=  betahat.ash.zero,
              betahat.ash.true=betahat.ash.true,
              pi0=pi0))
}

nonzerosim_unif=function(mixsd,mixpi_alt,betamean=5,bsd=1,method="shrink",mixcompdist="normal",minpi0=0,seedval=2014,nsamp=1000,niter=20){
  set.seed(seedval)
  beta =list()
  betahatsd=list()
  betahat = list()
  betahat.ash.nonzero = list()
  betahat.ash.zero = list()
  betahat.ash.true=list() 
  pi0 = rep(0,niter)	
  for(i in 1:niter){
    pi0[i]=runif(1,minpi0,1)
    mixpi = c(pi0[i],(1-pi0[i])*mixpi_alt)
    sds = sample(mixsd,nsamp,prob=mixpi,replace=TRUE)
    beta[[i]] = betamean+ sds*(2*runif(nsamp)-1)
    betahatsd[[i]] = bsd
    betahat[[i]] = beta[[i]]+rnorm(nsamp,0,betahatsd[[i]])
    betahat.ash.zero[[i]] = ash(betahat[[i]],betahatsd[[i]],method=method,mixcompdist=mixcompdist,nonzeromean=FALSE)
    betahat.ash.nonzero[[i]] = ash(betahat[[i]],betahatsd[[i]],method=method,mixcompdist=mixcompdist,nonzeromean=TRUE)
    g=unimix(mixpi,-mixsd,mixsd)
    betahat.ash.true[[i]] = ash(betahat[[i]]-betamean,betahatsd[[i]],method=method, mixcompdist="uniform",g=g)
    betahat.ash.true[[i]]$PosteriorMean=betahat.ash.true[[i]]$PosteriorMean+betamean
	if(class(betahat.ash.true[[i]]$fitted.g)=="normalmix"){
		betahat.ash.true[[i]]$fitted.g$mean=rep(betamean,length(betahat.ash.true[[i]]$fitted.g$pi))
	}else{
		betahat.ash.true[[i]]$fitted.g$a=betahat.ash.true[[i]]$fitted.g$a+betamean
		betahat.ash.true[[i]]$fitted.g$b=betahat.ash.true[[i]]$fitted.g$b+betamean
	}    
    print(betahat.ash.true[[i]]$PosteriorMean[1:20])
    print(i)
    toc()
}
  return(list(beta =beta,
              betahatsd=betahatsd,
              betahat = betahat,        
			  betahat.ash.nonzero=  betahat.ash.nonzero,
     		  betahat.ash.zero=  betahat.ash.zero,
              betahat.ash.true=betahat.ash.true,
              pi0=pi0))
}

nonzerosim_halfunif=function(mixsd,mixpi_alt,betamean=5,bsd=1,method="shrink",mixcompdist="normal",minpi0=0,seedval=2014,nsamp=1000,niter=20){
  set.seed(seedval)
  beta =list()
  betahatsd=list()
  betahat = list()
  betahat.ash.nonzero = list()
  betahat.ash.zero = list()
  betahat.ash.true=list() 
  pi0 = rep(0,niter)
  
  for(i in 1:niter){
    pi0[i]=runif(1,minpi0,1)
    mixpi = c(pi0[i],0.5*(1-pi0[i])*mixpi_alt,0.5*(1-pi0[i])*mixpi_alt)
    mixsdh=c(mixsd,-mixsd[-1])
    sds = sample(mixsdh,nsamp,prob=mixpi,replace=TRUE)
    beta[[i]] = betamean+ sds*(2*runif(nsamp)-1)
    betahatsd[[i]] = bsd
    betahat[[i]] = beta[[i]]+rnorm(nsamp,0,betahatsd[[i]])
    betahat.ash.zero[[i]] = ash(betahat[[i]],betahatsd[[i]],method=method,mixcompdist=mixcompdist,nonzeromean=FALSE)
    betahat.ash.nonzero[[i]] = ash(betahat[[i]],betahatsd[[i]],method=method,mixcompdist=mixcompdist,nonzeromean=TRUE)   
    g=unimix(mixpi,pmin(0,mixsdh),pmax(0,mixsdh))
    betahat.ash.true[[i]] = ash(betahat[[i]]-betamean,betahatsd[[i]],method=method, mixcompdist="halfuniform",g=g)
    betahat.ash.true[[i]]$PosteriorMean=betahat.ash.true[[i]]$PosteriorMean+betamean
	if(class(betahat.ash.true[[i]]$fitted.g)=="normalmix"){
		betahat.ash.true[[i]]$fitted.g$mean=rep(betamean,length(betahat.ash.true[[i]]$fitted.g$pi))
	}else{
		betahat.ash.true[[i]]$fitted.g$a=betahat.ash.true[[i]]$fitted.g$a+betamean
		betahat.ash.true[[i]]$fitted.g$b=betahat.ash.true[[i]]$fitted.g$b+betamean
	}    
    print(betahat.ash.true[[i]]$PosteriorMean[1:20])
    print(i)
    toc()
}
  return(list(beta =beta,
              betahatsd=betahatsd,
              betahat = betahat,        
			  betahat.ash.nonzero=  betahat.ash.nonzero,
     		  betahat.ash.zero=  betahat.ash.zero,
              betahat.ash.true=betahat.ash.true,
              pi0=pi0))
}


addingbackmean=function(betamean,sims){
	for(i in 1:length(betamean)){
	for(j in 1:length(sims[[i]]$betahat.ash.true)){
		sims[[i]]$betahat.ash.true[[j]]$PosteriorMean=sims[[i]]$betahat.ash.true[[j]]$PosteriorMean+betamean[i]
		if(class(sims[[i]]$betahat.ash.true[[j]]$fitted.g)=="normalmix"){
			sims[[i]]$betahat.ash.true[[j]]$fitted.g$mean=rep(betamean[i],length(sims[[i]]$betahat.ash.true[[j]]$fitted.g$pi))
		}else{
			sims[[i]]$betahat.ash.true[[j]]$fitted.g$a=sims[[i]]$betahat.ash.true[[j]]$fitted.g$a+betamean[i]
			sims[[i]]$betahat.ash.true[[j]]$fitted.g$b=sims[[i]]$betahat.ash.true[[j]]$fitted.g$b+betamean[i]
		}
	}
	}
	return(sims)
}


rmse = function(x,y){sqrt(mean((x-y)^2))} 
get_rmse.ash = function(a,b){rmse(a$PosteriorMean,b)}
plot_rmse_boxplot_nonzero = function(sims,scenarioid,inczero=FALSE,incbetahat=FALSE){
  res=list()
  for(i in 1:length(sims)){
    err.bayes = mapply(get_rmse.ash,sims[[i]]$betahat.ash.true,sims[[i]]$beta)
    err.ash.nonzero= mapply(get_rmse.ash,sims[[i]]$betahat.ash.nonzero,sims[[i]]$beta)
    err.ash.zero = mapply(get_rmse.ash,sims[[i]]$betahat.ash.zero,sims[[i]]$beta)
    err.betahat = mapply(rmse,sims[[i]]$betahat,sims[[i]]$beta)
    err.zero = unlist(lapply(sims[[i]]$beta,rmse,y=0))        
    res[[i]] = data.frame(Scenario= scenarioid[i],
                          ash.nonzero=err.ash.nonzero/err.bayes,
                          ash.zero= err.ash.zero/err.bayes
    )
    if(inczero){
      res[[i]]=data.frame(res[[i]],zero=err.zero/err.bayes)
    }
    if(incbetahat){
      res[[i]]=data.frame(res[[i]],betahat=err.betahat/err.bayes)
    }
  }
  res.melt = melt(res, id.vars=c("Scenario"),variable.name="Method")    
  ggplot(res.melt,aes(Method,value,color=Method)) + coord_flip() + geom_boxplot() + facet_grid(.~Scenario)
}



plot_loglik_boxplot_nonzero= function(sims, scenarioid){
  res=list()
  for(i in 1:length(sims)){
    loglik.bayes = mapply(get_loglik,sims[[i]]$betahat.ash.true)
    loglik.ash.nonzero= mapply(get_loglik,sims[[i]]$betahat.ash.nonzero)
    loglik.ash.zero = mapply(get_loglik,sims[[i]]$betahat.ash.zero)
    res[[i]] = data.frame(Scenario= scenarioid[i],
                          ash.nonzero= loglik.ash.nonzero-loglik.bayes,
                          ash.zero= loglik.ash.zero-loglik.bayes)   
  }
  require(reshape2)
  res.melt = melt(res, id.vars=c("Scenario"),variable.name="Method")    
  ggplot(res.melt,aes(Method,value,color=Method)) + geom_boxplot() + facet_grid(.~Scenario)
}


##Actual simulation
##on slave1
tic()
niter=100
nsamp=200
sim1=list()
sim2=list()
sim3=list()
betamean=c(0,0.1,2,10)
for(i in 1:length(betamean)){
	sim1[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="normal",niter= niter,nsamp= nsamp)
	sim2[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="uniform",niter= niter,nsamp= nsamp)
	sim3[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="halfuniform",niter= niter,nsamp= nsamp)
	cat("Iteration i=",i)
	toc()
}
save.image(file="simn123temp.RData")

##on slave2
tic()
niter=100
nsamp=200
sim4=list()
sim5=list()
sim6=list()
betamean=c(0,0.1,2,10)
for(i in 1:length(betamean)){
	sim4[[i]]=nonzerosim_unif(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="normal",niter= niter,nsamp= nsamp)
	sim5[[i]]=nonzerosim_unif(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="uniform",niter= niter,nsamp= nsamp)
	sim6[[i]]=nonzerosim_unif(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="halfuniform",niter= niter,nsamp=nsamp)
	cat("Iteration i=",i)
	toc()
}
save.image(file="simn456temp.RData")

##on slave3
tic()
niter=100
nsamp=200
sim7=list()
sim8=list()
sim9=list()
betamean=c(0,0.1,2,10)
for(i in 1:length(betamean)){
	sim7[[i]]=nonzerosim_halfunif(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="normal",niter=niter,nsamp=nsamp)
	sim8[[i]]=nonzerosim_halfunif(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="uniform",niter=niter,nsamp=nsamp)
	sim9[[i]]=nonzerosim_halfunif(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="halfuniform",niter=niter,nsamp=nsamp)
	cat("Iteration i=",i)
	toc()
}
save.image(file="simn789temp.RData")


##on slave4
tic()
set.seed(925)
niter=100
nsamp=1000
sim10=list()
betamean=c(0,0.5,1,2,3)
for(i in 1:length(betamean)){
	sim10[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],method="shrink",mixcompdist="normal",niter=niter,nsamp=nsamp)
	cat("Iteration i=",i)
	toc()
}
save.image(file="simn10temp.RData")



##Loading back the data and generate plots
load("simn123temp.RData")
load("simn456temp.RData")
load("simn789temp.RData")
load("simn10temp.RData")



pdf("sim1_normal.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim1,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Normal mixture")
dev.off()
pdf("sim2_uniform.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim2,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Uniform mixture")
dev.off()
pdf("sim3_halfuniform.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim3,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Halfuniform mixture")
dev.off()

pdf("sim4_normal.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim4,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Normal mixture")
dev.off()
pdf("sim5_uniform.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim5,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Uniform mixture")
dev.off()
pdf("sim6_halfuniform.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim6,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Halfuniform mixture")
dev.off()

pdf("sim7_normal.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim7,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Normal mixture")
dev.off()
pdf("sim8_uniform.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim8,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Uniform mixture")
dev.off()
pdf("sim9_halfuniform.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim9,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Halfuniform mixture")
dev.off()

pdf("sim10_normal.pdf",width=10,height=4)
plot_rmse_boxplot_nonzero(sim10,c("mean=0","mean=0.5","mean=1","mean=2","mean=3"))+ggtitle("Normal mixture")
dev.off()



pdf("sim1_normal_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim1,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Normal mixture")
dev.off()
pdf("sim2_uniform_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim2,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Uniform mixture")
dev.off()
pdf("sim3_halfuniform_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim3,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Halfuniform mixture")
dev.off()

pdf("sim4_normal_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim4,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Normal mixture")
dev.off()
pdf("sim5_uniform_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim5,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Uniform mixture")
dev.off()
pdf("sim6_halfuniform_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim6,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Halfuniform mixture")
dev.off()

pdf("sim7_normal_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim7,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Normal mixture")
dev.off()
pdf("sim8_uniform_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim8,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Uniform mixture")
dev.off()
pdf("sim9_halfuniform_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim9,c("mean=0","mean=0.1","mean=2","mean=10"))+ggtitle("Halfuniform mixture")
dev.off()

pdf("sim10_normal_loglik.pdf",width=10,height=4)
plot_loglik_boxplot_nonzero(sim10,c("mean=0","mean=0.5","mean=1","mean=2","mean=3"))+ggtitle("Normal mixture")
dev.off()

