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
  pvalue = 2*pnorm(-abs(zscore),lower.tail=TRUE)
  res = qvalue(pvalue)
  return(res)
}

qvalue2pi0_est = function(output){
  if(!is.list(output)){return(list(pi0_est=NA))} #deal with case where ERROR thrown
  else{
    return (list(pi0_est=output$pi0))
  }
}
