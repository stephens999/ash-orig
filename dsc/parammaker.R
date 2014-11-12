#' @title parammaker for FDR DSC
#'
#' @description Creates params for a DSC for methods to compute FDR/q values  
#' @details None
#' 
#' @param seed (integer) the seed for the pseudo-rng to be set before generating the parameters
#' #'
#' @return a list with the following elements
#' \item{nsamp}{Number of samples to use}
#' \item{g}{An object of class normalmix}
#' \item{dataseed}{the seed to use when generating the data}
#' \item{seed}{the seed used when generating the parameters ()}
#'
parammaker = function(seed){
  set.seed(seed)
  param = list(nsamp=1000,g=normalmix(c(0.5,0.5),c(0,0),c(0,1)), dataseed=seed+1,seed=seed)
  save(param,file=paramfile(seed))
  return(param)
}
