#' @title wrapper for ash for FDR DSC
#'
#' @description Runs ash to compute FDR/q values  
#' @details None
#' 
#' @param input a list with elements betahat and sebetahat
#' 
#' @return output a list with the following elements
#' \item{qvalue}{vector of qvalues, with jth element being the q value corresponding to (betahat_j,sebetahat_j)
#'
library(ashr)

ash.wrapper=function(input){
  res = ash(input$betahat,input$sebetahat,mixcompdist="halfuniform",method="fdr")
  return(list(qvalue=res$qvalue))
}