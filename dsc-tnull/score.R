#' @title compute score for shrinkage DSC
#'
#' @description Outputs the RMSE and MAE of estimated beta values 
#' @details None
#' 
#' @param data
#' @param output
#' 
#' @return score a list with
#' \item{RMSE}{root mean squared error of estimated beta values}
#' \item{MAE}{Median absolute error of estimated beta values}
#' \item{dFSR}{The number of discoveries at (tail) FSR=0.05}
#' \item{dFDR}{The number of discoveries at tail FDR = 0.05}
#'
score = function(data, output){
  return(list(RMSE=sqrt(mean((data$meta$beta-output$beta_est)^2)),
              MAE = median(abs(data$meta$beta-output$beta_est)),
              dFSR = sum(qval.from.lfdr(as.vector(output$lfsr))<0.05),
              dFDR = sum(output$q<0.05),
              pi0 = data$meta$pi0,
              pi0_est = output$pi0_est)
  )
}
