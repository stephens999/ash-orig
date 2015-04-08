null_t_datamaker = function(args){
  
   #assuming true standard deviation is 1, then betahatsd is rchisq
  sebetahat = sqrt(rchisq(args$nsamp, args$df)/args$df)
  
  beta = rep(0,args$nsamp)  
  betahat = beta + rnorm(args$nsamp,0,1)
  
  
  meta=list(beta=beta,pi0=1)
  input=list(betahat=betahat,sebetahat=sebetahat,df=args$df)
  
  #end of meat of function
  
  data = list(meta=meta,input=input)
  
  return(data)
  
}