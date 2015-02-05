#' @title wrapper for ash for shrinkage DSC
#'
#' @description Runs ash to compute betahat values  
#' @details None
#' 
#' @param input a list with elements betahat and sebetahat
#' @param args a list containing other additional arguments to ash
#' 
#' @return output a list with the following elements
#' \item{beta_est}{vector containing point estimates for beta}
#'
library(ashr)

ash.wrapper=function(input,args=NULL){
  if(is.null(args)){
    args=list(mixcompdist="halfuniform",method="fdr")
  }
  res = do.call(ash, args= c(list(betahat=input$betahat,sebetahat=input$sebetahat),args))
  return(list(beta_est=res$PosteriorMean,pi0_est=get_pi0(res)))
}
