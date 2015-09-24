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
  res = do.call(ash, args= c(list(betahat=input$betahat,sebetahat=input$sebetahat,df=input$df),args))
  return(res)
}

#uses ash function to compute bayes rule by running it with true g
bayes.wrapper = function(input, meta, args=NULL){
  pi0 = meta$pi0
  g1 = meta$g1 
  g=g1
  #create g by adding null component to g1
  g$pi= c(pi0, (1-pi0)*g1$pi)
  g$mean = c(0,g1$mean)
  g$sd = c(0,g1$sd)
  #do computations for bayes rule by running ash with g, with 0 iterations
  res = do.call(ash, args= list(betahat=input$betahat,sebetahat=input$sebetahat,df=input$df,g=g,control=list(maxiter=0)))
  return(res)
}

ash2beta_est =function(output){
  return (list(beta_est=output$PosteriorMean))
} 

ash2pi0_est =function(output){
  return (list(pi0_est=get_pi0(output)))
} 

ash2fitted.g = function(output){
  return (list(fitted.g=output$fitted.g))
}
