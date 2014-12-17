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
  zscore = input$betahat/input$sebetahat
  pvalue = pchisq(zscore^2,df=1,lower.tail=F)
  res = qvalue(pvalue)
  return(list(qvalue=res$qvalue,pi0=res$pi0))  
}