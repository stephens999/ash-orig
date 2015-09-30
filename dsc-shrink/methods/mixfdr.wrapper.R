#' @title wrapper for mixfdr for shrinkage DSC
#'
#' @description Runs mixfdr to compute effect estimates from z scores (computed from betahat and sebetahat)  
#' @details None
#' 
#' @param input a list with elements betahat and sebetahat
#' @param args a list containing other additional arguments to ash
#' 
#' @return output a list with the following elements
#' \item{beta_est}{vector containing point estimates for beta}
#' \item{pi0_est}{scalar estimate of pi0}
#'
library(mixfdr)

mixfdr.wrapper=function(input,args=NULL){
  res = try(do.call(mixFdr, args= c(list(x=input$betahat/input$sebetahat, noiseSD=1, plots=FALSE),args)))
  if(inherits(res,"try-error")){res=list(effectSize=rep(NA,length(input$betahat)),pi0_est=NA,pi=NA,mu=NA,sigma=NA,noiseSD=NA)}
  return(list(res=res,input=input))
}

mixfdr2beta_est = function(output){
  return (list(beta_est=output$res$effectSize*output$input$sebetahat))
}

mixfdr2pi0_est = function(output){
  return (list(pi0_est=output$res$pi[1]))
}

mixfdr2fitted.g = function(output){
  return (list(fitted.g=
    normalmix(output$res$pi,output$res$mu,
              sqrt(output$res$sigma^2-output$res$noiseSD^2))))
}

