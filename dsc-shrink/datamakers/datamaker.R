library(ashr)

#' @title datamaker for FDR DSC
#'
#' @description Simulates data for a DSC for methods to compute FDR/q values  
#' @details None
#' 
#' @param args A list of the remaining arguments, which in this case is
#' \item{nsamp}{The number of samples to create}
#' \item{g}{An object of class normalmix specifying the mixture distribution from which non-null beta values 
#' are to be simulated}
#' \item{min_pi0}{The minimum value of pi0, the proportion of true nulls}
#' \item{max_pi0}{The maximum value of pi0, the proportion of true null}
#' \item{betahatsd}{The standard deviation of betahat to use}
#'
#' @return a list with the following elements
#' \item{meta}{A list containing the meta data. In this case beta}
#' \item{input}{A list containing the input data; in this case the set of betahat values and their standard errors}
#' 
rnormmix_datamaker = function(args){
  #here is the meat of the function that needs to be defined for each dsc to be done
  pi0 = runif(1,args$min_pi0,args$max_pi0) #generate the proportion of true nulls randomly
  
  k = ncomp(args$g)
  comp = sample(1:k,args$nsamp,mixprop(args$g),replace=TRUE) #randomly draw a component
  isnull = (runif(args$nsamp,0,1) < pi0)
  beta = ifelse(isnull, 0,rnorm(args$nsamp,comp_mean(args$g)[comp],comp_sd(args$g)[comp]))
  sebetahat = args$betahatsd
  betahat = beta + rnorm(args$nsamp,0,sebetahat)
  meta=list(g1=args$g,beta=beta,pi0=pi0)
  input=list(betahat=betahat,sebetahat=sebetahat,df=NULL)
  
  #end of meat of function
  
  data = list(meta=meta,input=input)
  
  return(data)
  
}


