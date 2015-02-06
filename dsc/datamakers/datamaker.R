library(ashr)

#' @title datamaker for FDR DSC
#'
#' @description Simulates data for a DSC for methods to compute FDR/q values  
#' @details None
#' 
#' @param seed The seed for the pseudo-rng set before generating the parameters
#' @param args A list of the remaining arguments, which in this case is
#' \item{nsamp}{The number of samples to create}
#' \item{g}{An object of class normalmix specifying the mixture distribution from which beta values 
#' are to be simulated}
#'
#' @return a list with the following elements
#' \item{meta}{A list containing the meta data. In this case beta}
#' \item{input}{A list containing the input data for fdr methods; in this case the set of betahat values and their standard errors}
#' 
datamaker = function(seed,args){
  
  set.seed(seed)
  #here is the meat of the function that needs to be defined for each dsc to be done
  
  k = ncomp(args$g)
  comp = sample(1:k,args$nsamp,mixprop(args$g),replace=TRUE) #randomly draw a component
  beta = rnorm(args$nsamp,comp_mean(args$g)[comp],comp_sd(args$g)[comp])
  sebetahat = 1  
  betahat = beta + rnorm(args$nsamp,0,sebetahat)
  meta=list(beta=beta)
  input=list(betahat=betahat,sebetahat=sebetahat)
  
  #end of meat of function
  
  data = list(meta=meta,input=input)
  
  return(data)
  
}


