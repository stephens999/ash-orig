#' @title parammaker for FDR DSC
#'
#' @description Creates params for a DSC for methods to compute FDR/q values  
#' @details None
#' 
#' @param index a list with elements
#' \item{seed}{The seed for the pseudo-rng to be set before generating the parameters (may not affect parameters if no pseudo-rng used, but set anyway)}
#' \item{scenario}{the scenario to be used when generating parameters}
#' 
#' @return a list with the following elements
#' \item{nsamp}{Number of samples to use}
#' \item{g}{An object of class normalmix}
#' \item{dataseed}{the seed to use when generating the data}
#' \item{seed}{the seed used in the call} 
#'
parammaker = function(indexlist){
  set.seed(indexlist$seed)
  
  param = list(nsamp=1000,g=normalmix(c(0.5,0.5),c(0,0),c(0,3)), dataseed=indexlist$seed+1,seed=indexlist$seed)
  save(param,file=paramfile(indexlist))
  return(param)
}
