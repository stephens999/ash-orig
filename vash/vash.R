# Prior of se: se|pi,alpha.vec,c ~ pi*IG(alpha_k,c*(alpha_k-1))
# Likelihood: varhat|se ~ sj*Gamma(n/2,n/2)
# pi, alpha.vec, c: known
# Posterior weight of P(se|varhat) (IG mixture distn)
post_pi = function(n,varhat,alpha.vec,modalpha.vec,c,pi){
  N = length(varhat)
  K = length(alpha.vec)
  post.pi.mat = exp(outer(rep(1,N),log(pi))+n/2*log(n/2)-lgamma(n/2)
                  +(n/2-1)*outer(log(varhat),rep(1,K))
                  +outer(rep(1,N),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+n/2))
                  -outer(rep(1,N),alpha.vec+n/2)*log(outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  return(pimat=post.pi.mat)
}

# Normalize pi to make sum(pi)=1
normalized_pi = function(pi.mat){
  n = dim(pi.mat)[2]
  pi.normalized = pi.mat/outer(rowSums(pi.mat),rep(1,n))  
  return(pi.normalized)
}

# Posterior distn P(se|varhat)
# alpha.vec, c, pi: known
# varhat: observed (standard errors)^2
post_distn = function(n,varhat,alpha.vec,modalpha.vec,c,pi){ 
  N = length(varhat)
  K = length(alpha.vec)
  post.pi = normalized_pi(post_pi(n,varhat,alpha.vec,modalpha.vec,c,pi))
  
  post.IG.parama = outer(rep(1,N),alpha.vec+n/2)
  post.IG.paramb = outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))
  
  # Posterior mean: E(var|varhat)
  post.mixmean = post.IG.paramb/(post.IG.parama-1)
  post.mean = apply(post.pi*post.mixmean,1,sum)
  
  return(list(pi=post.pi,mean=post.mean,gammaa=post.IG.parama,gammab=post.IG.paramb))
}

# Log-likelihood: L(varhat|c,pi,alpha.vec)
loglike = function(N,K,logc,alpha.vec,n,varhat,pi,unimodal){
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  
  c=exp(logc)
  pimat = exp(outer(rep(1,N),log(pi))+n/2*log(n/2)-lgamma(n/2)
                  +(n/2-1)*outer(log(varhat),rep(1,K))
                  +outer(rep(1,N),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+n/2))
                  -outer(rep(1,N),alpha.vec+n/2)*log(outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  #classprob=pimat/rowSums(pimat)
  logl = sum(log(rowSums(pimat)))
  return(-logl)
}

# Gradient of funtion loglike (w.r.t logc)
gradloglik = function(N,K,logc,alpha.vec,n,varhat,pi,unimodal){
  c=exp(logc)
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  pimat = exp(outer(rep(1,N),log(pi))+n/2*log(n/2)-lgamma(n/2)
            +(n/2-1)*outer(log(varhat),rep(1,K))
            +outer(rep(1,N),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+n/2))
            -outer(rep(1,N),alpha.vec+n/2)*log(outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  classprob = pimat/rowSums(pimat)
  gradmat = c*classprob*(outer(rep(1,N),alpha.vec/c)
                     -outer(rep(1,N),(alpha.vec+n/2)*modalpha.vec)/
    (outer(rep(1,N),c*modalpha.vec)+outer(n/2*varhat,rep(1,K))))
  grad = sum(-gradmat)
  return(grad)
}

# Gradient of funtion loglike for single component prior (w.r.t logalpha)
gradloglik.a = function(N,K,logc,logalpha.vec,n,varhat,pi,unimodal){
  alpha.vec=exp(logalpha.vec)
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  c=exp(logc)
  grad=-alpha.vec*sum(logc+log(modalpha.vec)+alpha.vec/modalpha.vec-digamma(alpha.vec)+digamma(alpha.vec+n/2)
           -c*(alpha.vec+n/2)/(c*modalpha.vec+n/2*varhat)-log(c*modalpha.vec+n/2*varhat))
  return(grad)
}


# EM algorithm to estimate pi, c, and alpha(only when singlecomp==TRUE)
# prior: nullbiased: add weight to the null component at each iteration; (NEED TESTING!)
#        uniform: uniform weight
# ltol: tolerance of convergence
# maxiter: max number of iterations
EMest_pi = function(varhat,n,c,alpha.vec,pi,prior='uniform', ltol=0.0001, maxiter=5000, unimodal, SGD, singlecomp){
  K = length(pi)
  N = length(varhat)
  nullcomp = which.max(alpha.vec)
  if(is.null(c)){
    est.c = TRUE
    c=mean(varhat)
  }else{
    est.c = FALSE
  }
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
    
  if(prior=="nullbiased"){ 
    prior = rep(1,K)
    prior[nullcomp] = 10
  }else if(prior=="uniform"){
    prior = rep(1,K)
  }
  loglik = rep(NA,maxiter)
  
  mm = post_pi(n,varhat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  loglik[1] = sum(log(m.rowsum))
  classprob = mm/m.rowsum
  for(i in 2:maxiter){    
    if (est.c==TRUE & SGD==FALSE){
      est=nlminb(logc,loglike,gradloglik,N=N,K=K,alpha.vec=alpha.vec,n=n,varhat=varhat,pi=pi,unimodal=unimodal)
      c=exp(est$par[1])  
    }else if (est.c==TRUE & SGD==TRUE){
      logc=log(c)-0.0001/sqrt(i)*gradloglik(N,K,log(c),alpha.vec,n,varhat,pi,unimodal)
      c=exp(logc)
    } 
    if (singlecomp==TRUE){
      logalpha.vec=max(log(2),log(alpha.vec)-0.0001/(i^0.25)*gradloglik.a(N,K,log(c),log(alpha.vec),n,varhat,pi,unimodal))
      alpha.vec=exp(logalpha.vec)
      if(unimodal=='variance'){
        modalpha.vec=alpha.vec+1
      }else if(unimodal=='precision'){
        modalpha.vec=alpha.vec-1
      }
    }
    pi = colSums(classprob)+prior-1
    #pi=colSums(classprob)
    pi = ifelse(pi<0,0,pi); pi=pi/sum(pi);
    mm = post_pi(n,varhat,alpha.vec,modalpha.vec,c,pi)
    m.rowsum = rowSums(mm)
    loglik[i] = sum(log(m.rowsum))
    classprob = mm/m.rowsum
    if(abs(loglik[i]-loglik[i-1])<ltol) break;      
  }
  converged = (i< maxiter)
  niter = min(c(i,maxiter))
  return(list(pi=pi,classprob=classprob,loglik.final=loglik,converged=converged,niter=niter,c=c,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec))  
}

# Adaptive shrinkage for variances
# Unimodal: 'variance': variances~Mix IG with common mode c
#           'precision': precisions~Mix Gamma with common mode c
# SGD: use stochastic gradient descent to est hyperparams
# singlecomp: fit single component prior (est both prior mean and var params by EB)
vash = function(varhat,n,prior='uniform',unimodal='precision',alpha.vec=NULL,c=NULL,pi=NULL,SGD=TRUE,singlecomp=FALSE){
  if(is.null(alpha.vec)){
    alpha.vec = 2+2^seq(-3,10)
  }
  if(singlecomp==TRUE){
    alpha.vec=max(1/(mean(varhat)^2*sd(varhat)^2),2)
  }
  if(is.null(pi)){
    pi = rep(1,length(alpha.vec))/length(alpha.vec)
  }
    
  pifit = EMest_pi(varhat,n,c,alpha.vec,pi,prior, ltol=0.001, maxiter=5000, unimodal,SGD, singlecomp)
  post = post_distn(n,varhat,pifit$alpha.vec,pifit$modalpha.vec,pifit$c,pifit$pi)
  return(list(PosteriorMean=post$mean,pifit=pifit,post=post,alpha=alpha.vec,c=pifit$c,pi=pifit$pi,unimodal=unimodal))
}

# function to plot the Empirical Bayes prior
# xmax: plot density on (0,xmax)
vashEBprior=function(vashobj,xmax){
  xgrid=seq(0.0001,xmax,by=0.01)
  if(vashobj$unimodal=='variance'){
    EBprior.var.sep=dgamma(outer(1/xgrid,rep(1,length(vashobj$alpha))),
                           shape=outer(rep(1,length(xgrid)),vashobj$alpha),
                           rate=vashobj$c*outer(rep(1,length(xgrid)),vashobj$alpha+1))*outer(1/xgrid^2,rep(1,length(vashobj$alpha)))
    EBprior.var=rowSums(outer(rep(1,length(xgrid)),vashobj$pi)*EBprior.var.sep)
    EBprior.prec.sep=dgamma(outer(xgrid,rep(1,length(vashobj$alpha))),
                            shape=outer(rep(1,length(xgrid)),vashobj$alpha),
                            rate=vashobj$c*outer(rep(1,length(xgrid)),vashobj$alpha+1))                    
    EBprior.prec=rowSums(outer(rep(1,length(xgrid)),vashobj$pi)*EBprior.prec.sep)
  }else if (vashobj$unimodal=='precision'){
    EBprior.var.sep=dgamma(outer(1/xgrid,rep(1,length(vashobj$alpha))),
                           shape=outer(rep(1,length(xgrid)),vashobj$alpha),
                           rate=vashobj$c*outer(rep(1,length(xgrid)),vashobj$alpha-1))*outer(1/xgrid^2,rep(1,length(vashobj$alpha)))
    EBprior.var=rowSums(outer(rep(1,length(xgrid)),vashobj$pi)*EBprior.var.sep)
    EBprior.prec.sep=dgamma(outer(xgrid,rep(1,length(vashobj$alpha))),
                            shape=outer(rep(1,length(xgrid)),vashobj$alpha),
                            rate=vashobj$c*outer(rep(1,length(xgrid)),vashobj$alpha-1))                   
    EBprior.prec=rowSums(outer(rep(1,length(xgrid)),vashobj$pi)*EBprior.prec.sep)
  }
  return(list(xgrid=xgrid,EBprior.var=EBprior.var,EBprior.prec=EBprior.prec, 
              EBprior.var.sep=EBprior.var.sep, EBprior.prec.sep=EBprior.prec.sep))
}