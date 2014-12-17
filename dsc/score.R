#' @title compute score for FDR DSC
#'
#' @description Outputs the actual FDR and positive FDR for simulated data at threshold 0.05  
#' @details None
#' 
#' @param param
#' @param data
#' @param output
#' 
#' @return score a list with
#' \item{fd}{number of false discoveries}
#' \item{td}{number of true discoveries}
#'
#'
score = function(data, output){
  return(list(fd=sum(output$qvalue[data$meta$beta==0]<0.05),
              td = sum(output$qvalue[data$meta$beta!=0]<0.05),
              pi0 = output$pi0))
}