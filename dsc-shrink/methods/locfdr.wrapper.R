#' @title wrapper for locfdr 
#'
#' @description Runs locfdr on z scores
#' @details None
#' 
#' @param input a list with elements betahat and sebetahat
#' @param args a list containing other additional arguments to locfdr
#' 
#' @return output a list 
#'
library(locfdr)

locfdr.wrapper=function(input,args=NULL){
  res = try(do.call(locfdr, args= c(list(zz=input$betahat/input$sebetahat, nulltype=0, plot=0),args)))
  if(inherits(res,"try-error")){res=list(fdr=rep(NA,length(input$betahat)),fp0=matrix(NA,nrow=6,ncol=3))}
  return(res)
}


locfdr2pi0_est = function(output){
  return (list(pi0_est=output$fp0[1,3]))
}


