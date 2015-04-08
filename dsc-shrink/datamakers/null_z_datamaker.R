null_z_datamaker = function(args){
  
  #assuming true standard deviation is 1, then betahatsd is rchisq
  sebetahat = rep(1,args$nsamp)
  
  beta = rep(0,args$nsamp)  
  betahat = beta + rnorm(args$nsamp,0,1)
  
  
  meta=list(beta=beta,pi0=1)
  input=list(betahat=betahat,sebetahat=sebetahat,df=NULL)
  
  #end of meat of function
  
  data = list(meta=meta,input=input)
  
  return(data)
  
}