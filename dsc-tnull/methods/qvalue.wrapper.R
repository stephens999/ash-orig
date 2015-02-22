#' @title wrapper for qvalue for FDR DSC
#'
#' @description Runs qvalue to compute FDR/q values  
#' @details None
#' 
#' @param input a list with elements betahat and sebetahat
#' @param add.args a list with additional arguments to qvalue 
#' 
#' @return output a list with the following elements
#' \item{qvalue}{vector of qvalues, with jth element being the q value corresponding to (betahat_j,sebetahat_j)
#'
#'
library(qvalue)

qvalue.wrapper = function(input,args=NULL){
  tscore = input$betahat/input$sebetahat
  pvalue = 2*pt(-abs(tscore),df=input$df,lower.tail=TRUE)
  res = qvalue(pvalue)
  return(list(beta_est=rep(NA,length(pvalue)),pi0_est=res$pi0,q=res$q,lfsr=rep(NA,length(pvalue))))  
}