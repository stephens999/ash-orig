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
#'
#'
score = function(data, output){
  return(list(RMSE=sqrt(mean((data$meta$beta-output$beta_est)^2)),
              MAE = median(abs(data$meta$beta-output$beta_est)),
              pi0 = data$meta$pi0,
              pi0_est = output$pi0_est)
  )
}