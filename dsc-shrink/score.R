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
              MAE = median(abs(data$meta$beta-output$beta_est))))
}

score2 = function(data, output){
  return(list(pi0 = data$meta$pi0,
              pi0_est = output$pi0_est))
}

score3 = function(data, output){
  return(c(S=pcdf_post(output$fitted.g,data$meta$beta,
              data$input$betahat,data$input$sebetahat,v=NULL)))
}

score_neg = function(data, output){
  return(c(S=output$NegativeProb))
}

score_pos = function(data, output){
  return(c(S=output$PositiveProb))
}

score_fdr = function(data, output){
  return(c(S=output$res$fdr))
}

score_lfsr = function(data, output){
  return(c(S=output$lfsr))
}

score_lfdr = function(data, output){
  return(c(S=output$lfdr))
}

score_betahat = function(data, output){
  return(c(S=data$input$betahat))
}

score_logLR = function(data,output){
  return(c(logLR=output$logLR))
}
