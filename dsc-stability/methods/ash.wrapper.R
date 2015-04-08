#' @title wrapper for ash for shrinkage DSC
#'
#' @description Runs ash to compute betahat values  
#' @details None
#' 
#' @param input a list with elements betahat and sebetahat
#' @param args a list containing other additional arguments to ash
#' 
#' @return output a list containing a vector of loglikelihoods from 10 independent runs of ash
#'
library(ashr)

ash.wrapper=function(input,args=NULL){
  if(is.null(args)){
    args=list(mixcompdist="halfuniform")
  }
  loglik=rep(0,11)
  
  #first run is with non-random start (ie default start)
  res=do.call(ash, args= c(list(betahat=input$betahat,sebetahat=input$sebetahat,randomstart=FALSE),args))
  loglik[1]= get_loglik(res)
    
  for(i in 2:11){
    res = do.call(ash, args= c(list(betahat=input$betahat,sebetahat=input$sebetahat,randomstart=TRUE),args))
    loglik[i]= get_loglik(res)
  }
  return(list(loglik=loglik))
}
