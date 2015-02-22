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
  res = do.call(ash, args= c(list(betahat=input$betahat,sebetahat=input$sebetahat,df = input$df),args))
  return(list(beta_est=res$PosteriorMean,pi0_est=get_pi0(res),lfsr=res$lfsr,q=res$qvalue))
}

#this wrapper deals with df by adjusting betahat to make betahat/sebetahat to normal first
ash.betaeff.wrapper=function(input,args=NULL){
  if(is.null(args)){
    args=list(mixcompdist="halfuniform",method="fdr")
  }
  #first convert beta
  betahat_eff = ashr:::effective.effect(input$betahat,input$sebetahat,input$df)
  res = do.call(ash, args= c(list(betahat=betahat_eff,sebetahat=input$sebetahat),args))
  return(list(beta_est=res$PosteriorMean,pi0_est=get_pi0(res),lfsr=res$lfsr,q=res$qvalue))
}

#this wrapper deals with df by adjusting sebetahat to make betahat/sebetahat normal first
ash.sigmaeff.wrapper=function(input,args=NULL){
  if(is.null(args)){
    args=list(mixcompdist="halfuniform",method="fdr")
  }
  #first convert beta
  betahat_eff = ashr:::effective.effect(input$betahat,input$sebetahat,input$df)
  se_eff = input$sebetahat * input$betahat/betahat_eff
  res = do.call(ash, args= c(list(betahat=input$betahat,sebetahat=se_eff),args))
  return(list(beta_est=res$PosteriorMean,pi0_est=get_pi0(res),lfsr=res$lfsr,q=res$qvalue))
}