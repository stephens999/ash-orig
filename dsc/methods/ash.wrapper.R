#' @title wrapper for ash for FDR DSC
#'
#' @description Runs ash to compute FDR/q values  
#' @details None
#' 
#' @param input a list with elements betahat and sebetahat
#' @param add.args a list containing other additional arguments to ash
#' 
#' @return output a list with the following elements
#' \item{qvalue}{vector of qvalues, with jth element being the q value corresponding to (betahat_j,sebetahat_j)
#'
library(ashr)

ash.wrapper=function(input,args=NULL){
  if(is.null(args)){
    args=list(mixcompdist="halfuniform",method="fdr")
  }
  res = do.call(ash, args= c(list(betahat=input$betahat,sebetahat=input$sebetahat),
                             args))
  return(list(qvalue=res$qvalue,pi0=get_pi0(res)))
}
