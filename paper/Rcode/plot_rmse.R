require(reshape2)
require(ggplot2)

rmse = function(x,y){sqrt(mean((x-y)^2))}
get_rmse.ash = function(a,b){rmse(a$PosteriorMean,b)}
get_rmse.mixfdr = function(a,b){rmse(a$effectSize,b)}
get_loglik.mixfdr = function(a,betahat,betahatsd){loglik_conv(normalmix(a$pi,a$mu,a$sigma-1),betahat,betahatsd,NULL)}
cdf.mixfdr = function(a,x){mixcdf(normalmix(a$pi,a$mu,a$sigma-1),x)}

plot_loglik_boxplot = function(sims){
  res=list()
  for(i in 1:length(sims)){
    loglik.bayes = mapply(get_loglik,sims[[i]]$fit.ash.true)
    loglik.ash.n= mapply(get_loglik,sims[[i]]$fit.ash.n)
    loglik.ash.u = mapply(get_loglik,sims[[i]]$fit.ash.u)
    loglik.ash.hu = mapply(get_loglik,sims[[i]]$fit.ash.hu)
    
    loglik.ash.fdr.n = mapply(get_loglik,sims[[i]]$fit.ash.fdr.n)         
    res[[i]] = data.frame(Scenario=i,
                          ash.normal=loglik.ash.n-loglik.bayes,
                          ash.uniform=loglik.ash.u-loglik.bayes,
                          ash.halfuniform=loglik.ash.hu-loglik.bayes,
                          ash.normal.fdr = loglik.ash.fdr.n-loglik.bayes)
    
  }
  require(reshape2)
  res.melt = melt(res, id.vars=c("Scenario"),variable.name="Method")    
  ggplot(res.melt,aes(Method,value,color=Method)) + geom_boxplot() + facet_grid(.~Scenario)
  
}

plot_rmse = function(sims,inczero=FALSE,incbetahat=FALSE){
  res=list()
  
  for(i in 1:length(sims)){
    err.bayes = mapply(get_rmse.ash,sims[[i]]$fit.ash.true,sims[[i]]$beta)
    err.ash.n= mapply(get_rmse.ash,sims[[i]]$fit.ash.n,sims[[i]]$beta)
    err.ash.u = mapply(get_rmse.ash,sims[[i]]$fit.ash.u,sims[[i]]$beta)
    err.ash.hu = mapply(get_rmse.ash,sims[[i]]$fit.ash.hu,sims[[i]]$beta)
    err.mixfdr=mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr,sims[[i]]$beta)
    err.mixfdr.enull = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.enull,sims[[i]]$beta)
    err.mixfdr.J10 = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.J10,sims[[i]]$beta)
    #err.mixfdr.J100 = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.J100,sims[[i]]$beta)
    err.ash.fdr.n = mapply(get_rmse.ash,sims[[i]]$fit.ash.fdr.n,sims[[i]]$beta)
    
    err.betahat = mapply(rmse,sims[[i]]$betahat,sims[[i]]$beta)
    err.zero = unlist(lapply(sims[[i]]$beta,rmse,y=0))        
    res[[i]] = data.frame(Scenario=i,bayes=err.bayes,ash.normal=err.ash.n,ash.uniform=err.ash.u,ash.halfuniform=err.ash.hu,mixfdr=err.mixfdr,
                          mixfdr.enull = err.mixfdr.enull,
                          mixfdr.J10 = err.mixfdr.J10,
                          #mixfdr.J100 = err.mixfdr.J100,
                          ash.normal.fdr = err.ash.fdr.n)
    if(inczero){
      res[[i]]=data.frame(res[[i]],zero=err.zero)
    }
    if(incbetahat){
      res[[i]]=data.frame(res[[i]],betahat=err.betahat)
    }
  }

  res.melt = melt(res, id.vars=c("Scenario","bayes"),variable.name="Method")    
  p=ggplot(data=res.melt,aes(bayes,value,colour=Method)) +geom_point(shape=16) +
    facet_grid(. ~ Scenario,scale="free_x") +
    geom_abline(colour = "black") +
    xlab("RMSE (Optimal Bayes Rule)") +
    ylab("RMSE (other method)")
  print(p +
          coord_equal(ratio=1))
}




plot_rmse_boxplot = function(sims,inczero=FALSE,incbetahat=FALSE,incmixfdr=FALSE){
  res=list()
  
  for(i in 1:length(sims)){
    err.bayes = 1; #mapply(get_rmse.ash,sims[[i]]$fit.ash.true,sims[[i]]$beta)
    err.ash.n= mapply(get_rmse.ash,sims[[i]]$fit.ash.n,sims[[i]]$beta)
    err.ash.u = mapply(get_rmse.ash,sims[[i]]$fit.ash.u,sims[[i]]$beta)
    err.ash.hu = mapply(get_rmse.ash,sims[[i]]$fit.ash.hu,sims[[i]]$beta)
    err.mixfdr=mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr,sims[[i]]$beta)
    err.mixfdr.enull = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.enull,sims[[i]]$beta)
    err.mixfdr.J10 = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.J10,sims[[i]]$beta)
    #       err.mixfdr.J100 = mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr.J100,sims[[i]]$beta)
    err.ash.fdr.n = mapply(get_rmse.ash,sims[[i]]$fit.ash.fdr.n,sims[[i]]$beta)         
    err.betahat = mapply(rmse,sims[[i]]$betahat,sims[[i]]$beta)
    err.zero = unlist(lapply(sims[[i]]$beta,rmse,y=0))        
    res[[i]] = data.frame(Scenario=i,
                          ash.normal=err.ash.n/err.bayes,
                          ash.uniform=err.ash.u/err.bayes,
                          ash.halfuniform=err.ash.hu/err.bayes,
                          ash.normal.fdr = err.ash.fdr.n/err.bayes
    )
    if(incmixfdr){
      res[[i]]=data.frame(res[[i]],mixfdr=err.mixfdr/err.bayes,
                          mixfdr.enull = err.mixfdr.enull/err.bayes,
                          mixfdr.J10 = err.mixfdr.J10/err.bayes)
      #mixfdr.J100 = err.mixfdr.J100/err.bayes)
    }
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




plot_rmse_loglik_boxplot = function(sims){
  res=list()
  
  for(i in 1:length(sims)){
    err.bayes = mapply(get_rmse.ash,sims[[i]]$fit.ash.true,sims[[i]]$beta)
    err.ash.n= mapply(get_rmse.ash,sims[[i]]$fit.ash.n,sims[[i]]$beta)
    err.ash.u = mapply(get_rmse.ash,sims[[i]]$fit.ash.u,sims[[i]]$beta)
    err.ash.hu = mapply(get_rmse.ash,sims[[i]]$fit.ash.hu,sims[[i]]$beta)
    err.mixfdr=mapply(get_rmse.mixfdr,sims[[i]]$fit.mixfdr,sims[[i]]$beta)
    err.ash.fdr.n = mapply(get_rmse.ash,sims[[i]]$fit.ash.fdr.n,sims[[i]]$beta)         
    
    res[[i]] = data.frame(Scenario=i,
                          ash.normal=err.ash.n/err.bayes,
                          ash.uniform=err.ash.u/err.bayes,
                          ash.halfuniform=err.ash.hu/err.bayes,
                          ash.normal.fdr = err.ash.fdr.n/err.bayes)
    
  }
  
  res2=list()
  for(i in 1:length(sims)){
    loglik.bayes = mapply(get_loglik,sims[[i]]$fit.ash.true)
    loglik.ash.n= mapply(get_loglik,sims[[i]]$fit.ash.n)
    loglik.ash.u = mapply(get_loglik,sims[[i]]$fit.ash.u)
    loglik.ash.hu = mapply(get_loglik,sims[[i]]$fit.ash.hu)
    loglik.mixfdr=mapply(get_loglik.mixfdr,sims[[i]]$fit.mixfdr,sims[[i]]$betahat,sims[[i]]$betahatsd)
    loglik.ash.fdr.n = mapply(get_loglik,sims[[i]]$fit.ash.fdr.n)         
    res2[[i]] = data.frame(Scenario=i,
                           ash.normal=loglik.ash.n-loglik.bayes,
                           ash.uniform=loglik.ash.u-loglik.bayes,
                           ash.halfuniform=loglik.ash.hu-loglik.bayes,
                           ash.normal.fdr = loglik.ash.fdr.n-loglik.bayes)
    
  }
  res.melt = melt(res, id.vars=c("Scenario"),variable.name="Method")    
  res2.melt = melt(res2, id.vars=c("Scenario"),variable.name="Method")    
  res.melt$type='RMSE (vs Bayes Rule)'
  res2.melt$type='log(likelihood) (vs Bayes Rule)'
  ggplot(rbind(res.melt,res2.melt),aes(Method,value)) + geom_boxplot() + facet_grid(type~Scenario,scale="free_y")
  
}
