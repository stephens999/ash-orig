#' @title datamaker for FDR DSC
#'
#' @description Simulates data for a DSC for methods to compute FDR/q values  
#' @details None
#' 
#' @param indexlist a list with elements
#' \item{seed}{The seed for the pseudo-rng set before generating the parameters}
#' \item{scenario}{the scenario used when generating parameters}
#' This item points the datamaker to the correct parameter file - it indexes the parameter file
#' created by parammaker(indexlist)
#' 
#' \item{nsamp}{The number of samples to create}
#' \item{g}{An object of class normalmix specifying the mixture distribution from which beta values 
#' are to be simulated}
#'
#' @return a list with the following elements
#' \item{meta}{A list containing the meta data. In this case beta}
#' \item{input}{A list containing the input data for fdr methods; in this case the set of beta values and their standard errors}
#' 
datamaker = function(indexlist){
  load(file=paramfile(indexlist))
  set.seed(param$dataseed)
  
  k = ncomp(param$g)
  comp = sample(1:k,param$nsamp,mixprop(param$g),replace=TRUE) #randomly draw a component
  beta = rnorm(param$nsamp,comp_mean(param$g)[comp],comp_sd(param$g)[comp])
  sebetahat = 1  
  betahat = beta + rnorm(param$nsamp,0,sebetahat)
  meta=list(beta=beta)
  input=list(betahat=betahat,sebetahat=sebetahat)
  
  data = list(meta=meta,input=input)
  save(data,file=datafile(indexlist))
  return(data)
}