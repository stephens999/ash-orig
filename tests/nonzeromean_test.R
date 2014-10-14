# on Terminal
# open -n /Applications/RStudio.app
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

nonzerosim=function(mixsd,mixpi_alt,betamean=5,bsd=1,datacompdist="normal",datadf=NULL,method="shrink",mixcompdist="normal",minpi0=0,seedval=2014,nsamp=1000,niter=20){
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
    if(datacompdist=="normal"){
	  mixpi = c(pi0[i],(1-pi0[i])*mixpi_alt)
      sds = sample(mixsd,nsamp,prob=mixpi,replace=TRUE)
      beta[[i]] = rnorm(nsamp,betamean,sds)
    }else if(datacompdist=="uniform"){
      mixpi = c(pi0[i],(1-pi0[i])*mixpi_alt)
      sds = sample(mixsd,nsamp,prob=mixpi,replace=TRUE)
      beta[[i]] = betamean+ sds*(2*runif(nsamp)-1)   	
    }else if(datacompdist=="halfuniform"){
      mixpi = c(pi0[i],0.5*(1-pi0[i])*mixpi_alt,0.5*(1-pi0[i])*mixpi_alt)
      mixsdh=c(mixsd,-mixsd[-1])
      sds = sample(mixsdh,nsamp,prob=mixpi,replace=TRUE)
      beta[[i]] = betamean+ sds*runif(nsamp)
    }
    else{
    	stop("Invalid type of datacompdist supplied")
    }
    
    betahatsd[[i]] = bsd
    if(is.null(datadf)){
      betahat[[i]] = beta[[i]]+rnorm(nsamp,0,betahatsd[[i]])
    }else{
      betahat[[i]] = beta[[i]]+rt(nsamp,df=datadf)*betahatsd[[i]]
    }
    
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
    print(i)
  }
  toc()
  return(list(beta =beta,
              betahatsd=betahatsd,
              betahat = betahat,        
			  betahat.ash.nonzero=  betahat.ash.nonzero,
     		  betahat.ash.zero=  betahat.ash.zero,
              betahat.ash.true=betahat.ash.true,
              pi0=pi0))
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
    cc=!(err.bayes==0)
	a=log(err.ash.nonzero[cc]/err.bayes[cc])/log(2)
	b=log(err.ash.zero[cc]/err.bayes[cc])/log(2)
	   
    res[[i]] = data.frame(Scenario= scenarioid[i],
                          ash.nonzero=a,
                          ash.zero= b
    )
    if(inczero){
      res[[i]]=data.frame(res[[i]],zero=err.zero/err.bayes)
    }
    if(incbetahat){
      res[[i]]=data.frame(res[[i]],betahat=err.betahat/err.bayes)
    }
  }
  res.melt = melt(res, id.vars=c("Scenario"),variable.name="Method")    
  ggplot(res.melt,aes(Method,value,color=Method)) + coord_flip() + geom_boxplot() + facet_grid(.~Scenario)+labs(y="log2(RMSE/RMSE_BAYES)")
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
  ggplot(res.melt,aes(Method,value,color=Method)) + coord_flip() +geom_boxplot() + facet_grid(.~Scenario)+labs(y="loglikelihood(ash)-loglikelihood(Full Bayesian)")
}




##Actual simulation
##on slave1
tic()
niter=100
nsamp=200
sim1=list()
sim2=list()
sim3=list()
betamean=c(0,0.05,0.1,0.15,0.2)
for(i in 1:length(betamean)){
	sim1[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist="normal",method="shrink",mixcompdist="normal",niter= niter,nsamp= nsamp)
	sim2[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist="normal",method="shrink",mixcompdist="uniform",niter= niter,nsamp= nsamp)
	sim3[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist="normal",method="shrink",mixcompdist="halfuniform",niter= niter,nsamp= nsamp)
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
betamean=c(0,0.01,0.02,0.03,0.04)
datacompdist="uniform"
for(i in 1:length(betamean)){
	sim4[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist=datacompdist,method="shrink",mixcompdist="normal",niter= niter,nsamp= nsamp)
	sim5[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist=datacompdist,method="shrink",mixcompdist="uniform",niter= niter,nsamp= nsamp)
	sim6[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist=datacompdist,method="shrink",mixcompdist="halfuniform",niter= niter,nsamp= nsamp)
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
betamean=c(0,0.001,0.002,0.003,0.004)
datacompdist="uniform"
for(i in 1:length(betamean)){
	sim7[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist=datacompdist,method="shrink",mixcompdist="normal",niter= niter,nsamp= nsamp)
	sim8[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist=datacompdist,method="shrink",mixcompdist="uniform",niter= niter,nsamp= nsamp)
	sim9[[i]]=nonzerosim(mixsds,mixpi_alt,betamean[i],datacompdist=datacompdist,method="shrink",mixcompdist="halfuniform",niter= niter,nsamp= nsamp)
	cat("Iteration i=",i)
	toc()
}
save.image(file="simn789temp.RData")




##Loading back the data and generate plots
load("simn123temp.RData")
load("simn456temp.RData")
load("simn789temp.RData")
library(ggplot2)
library(reshape2)


scenarioid=c("mean=0","mean=0.05","mean=0.1","mean=0.15","mean=0.2")
pdf("sim1_normal.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim1, scenarioid)+ggtitle("Normal mixture data fitted using Normal mixture")
dev.off()
pdf("sim2_uniform.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim2,scenarioid)+ggtitle("Normal mixture data fitted using Uniform mixture")
dev.off()
pdf("sim3_halfuniform.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim3,scenarioid)+ggtitle("Normal mixture data fitted using Halfuniform mixture")
dev.off()

scenarioid=c("mean=0","mean=0.01","mean=0.02","mean=0.03","mean=0.04")
pdf("sim4_normal.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim4,scenarioid)+ggtitle("Uniform mixture data fitted using Normal mixture")
dev.off()
pdf("sim5_uniform.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim5,scenarioid)+ggtitle("Uniform mixture data fitted using Uniform mixture")
dev.off()
pdf("sim6_halfuniform.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim6,scenarioid)+ggtitle("Uniform mixture data fitted using Halfuniform mixture")
dev.off()

scenarioid=c("mean=0","mean=0.001","mean=0.002","mean=0.003","mean=0.004")
pdf("sim7_normal.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim7,scenarioid)+ggtitle("Halfuniform mixture data fitted using Normal mixture")
dev.off()
pdf("sim8_uniform.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim8,scenarioid)+ggtitle("Halfuniform mixture data fitted using Uniform mixture")
dev.off()
pdf("sim9_halfuniform.pdf",width=20,height=8)
plot_rmse_boxplot_nonzero(sim9,scenarioid)+ggtitle("Halfuniform mixture data fitted using Halfuniform mixture")
dev.off()

scenarioid=c("mean=0","mean=0.05","mean=0.1","mean=0.15","mean=0.2")
pdf("sim1_normal_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim1,scenarioid)+ggtitle("Normal mixture data fitted using Normal mixture")
dev.off()
pdf("sim2_uniform_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim2,scenarioid)+ggtitle("Normal mixture data fitted using Uniform mixture")
dev.off()
pdf("sim3_halfuniform_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim3,scenarioid)+ggtitle("Normal mixture data fitted using Halfuniform mixture")
dev.off()

scenarioid=c("mean=0","mean=0.01","mean=0.02","mean=0.03","mean=0.04")
pdf("sim4_normal_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim4,scenarioid)+ggtitle("Uniform mixture data fitted using Normal mixture")
dev.off()
pdf("sim5_uniform_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim5,scenarioid)+ggtitle("Uniform mixture data fitted using Uniform mixture")
dev.off()
pdf("sim6_halfuniform_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim6,scenarioid)+ggtitle("Uniform mixture data fitted using Halfuniform mixture")
dev.off()

scenarioid=c("mean=0","mean=0.001","mean=0.002","mean=0.003","mean=0.004")
pdf("sim7_normal_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim7,scenarioid)+ggtitle("Halfuniform mixture data fitted using Normal mixture")
dev.off()
pdf("sim8_uniform_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim8,scenarioid)+ggtitle("Halfuniform mixture data fitted using Uniform mixture")
dev.off()
pdf("sim9_halfuniform_loglik.pdf",width=20,height=8)
plot_loglik_boxplot_nonzero(sim9,scenarioid)+ggtitle("Halfuniform mixture data fitted using Halfuniform mixture")
dev.off()



sims=list(sim1,sim2,sim3,sim4,sim5,sim6,sim7,sim8,sim9)

rmsetable=function(sims, betamean){
  res=matrix(NA,nrow=length(sims),ncol=length(betamean))
  for(i in 1:length(sims)){
    err.bayes = mapply(get_rmse.ash,sims[[i]]$betahat.ash.true,sims[[i]]$beta)
    err.ash.nonzero= mapply(get_rmse.ash,sims[[i]]$betahat.ash.nonzero,sims[[i]]$beta)
    err.ash.zero = mapply(get_rmse.ash,sims[[i]]$betahat.ash.zero,sims[[i]]$beta)
    cc=!(err.bayes==0)
	a=log(err.ash.nonzero[cc]/err.bayes[cc])/log(2)
	b=log(err.ash.zero[cc]/err.bayes[cc])/log(2)
    res[i,1:length(betamean)] = c(betamean[i],mean(a),sd(a)/length(a),mean(b),sd(b)/length(b))
  }
  return(res)
}

##Generating Tables summarizing the result in markdown format
for(i in 1:9){
	if(i<4){
		betamean=c(0,0.05,0.1,0.15,0.2)
	}
	else if(i<7){
		betamean=c(0,0.01,0.02,0.03,0.04)
	}
	else{
		betamean=c(0,0.001,0.002,0.003,0.004)
	}
	a=rmsetable(sims[[i]],betamean)
	colnames(a)=c("True Mean","log2(RMSE_nonzero/RMSE_BAYES)","SE","log2(RMSE_zero/RMSE_BAYES)","SE")
	kable(a,format="markdown",digits=4,align='c')
}

