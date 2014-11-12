#' @title return the path to a data file, parameter file, output file or results file
#'
#' @description return the path to a data file, parameter file, output file or results file
#' 
#' @param seed (integer) the seed for the pseudo-rng which identifies/indexes the file
#' @param datadir/paramdir/outputdir/resultsdir the (relative) path to the directory containing the relevant files
#' 
#' @return string containing path to file
#' 
#' @example datafile(1)
#' 

datafile = function(seed,datadir="data"){
  return(file.path(datadir,paste0("data.",seed,".RData")))
}

paramfile = function(seed,paramdir="param"){
  return(file.path(paramdir,paste0("param.",seed,".RData")))
}

outputfile = function(methodname,seed,outputdir="output"){
  return(file.path(outputdir,methodname,paste0(methodname,".output.",seed,".RData")))
}

resultsfile = function(methodname,seed,resultsdir="results"){
  return(file.path(resultsdir,methodname,paste0(methodname,".results.",seed,".RData")))
}

#' @title Apply a method to an input and produce output
#'
#' @description Apply a single method to a single input trial and produce (and save) corresponding output
#' 
#' @param seed (integer) the seed for the pseudo-rng which identifies/indexes trials
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' 
#' @return output a list of appropriate format to be determined by the comparison being run
#' 
#' 
apply_method_singletrial=function(seed,method){
  load(file=datafile(seed))
  output=do.call(method$fn,list(input=data$input))
  save(output,file=outputfile(method$name,seed))
  return(output)
}

#' @title Apply a method to all inputs and produce corresponding outputs
#'
#' @description Apply a method to all inputs and produce corresponding outputs (by repeated application of apply_method_once).
#' Results are saved in a file
#'  
#' @param seed (list of integers) the seeds for the pseudo-rng which identifies/indexes trials
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' 
#' @return none - the results are saved in a file determined by outputfile(seed, method$methodname)
#' 
apply_method = function(seed,method){
  l_ply(seed,apply_method_singletrial,method=method)
}

#' @title Score a method on a single trial and save results
#'
#' @description Score results of a single method for a single trial and produce (and save) corresponding results
#' 
#' @param seed (integer) the seed for the pseudo-rng which identifies/indexes trials
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' @param scorefn a function that scores output based on comparisons with input, parameters and metadata
#'
#' @return results, a list of appropriate format to be determined by the comparison being run (maybe required to be a dataframe?)
#' 
#' 
score_method_singletrial = function(seed,method,scorefn){
  load(file=paramfile(seed))
  load(file=datafile(seed))
  load(file=outputfile(method$name,seed))
  
  results=scorefn(param,data,output)
  save(results,file=resultsfile(method$name,seed))
  return(results)
}

#' @title Score a method on all trials and save results
#'
#' @description Score a method on all trials and save results (by repeated application of score_method_singletrial)
#' 
#' @param seed (list of integers) the seeds for the pseudo-rng which identifies/indexes trials
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' @param scorefn a function that scores output based on comparisons with input, parameters and metadata
#'
#' @return results, a list of appropriate format to be determined by the comparison being run (maybe required to be a dataframe?)
#' 
#'
score_method=function(seed,method,scorefn){
  l_ply(seed,score_method_singletrial,method=method,scorefn=scorefn)
}





#' @title Get the results of a single method for a single trial
#'
#' @description Get the results of a single method for a single trial
#' 
#' @param seed (integer) the seed for the pseudo-rng which identifies/indexes trials
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#'
#' @return results, a data frame of results, the details will depend on the comparison being run
#' 
get_one_result = function(seed,method,resultsdir="results"){
  load(file=resultsfile(method$name,seed))
  return(data.frame(results))
}

#' @title Get the results of a single method for multiple trials
#'
#' @description Get the results of a single method for multiple trials
#' 
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' @param seed (list of integers) the seeds for the pseudo-rng which identifies/indexes trials
#' @return a data frame of results, with one row for each trial. The details of the columns will depend on the comparison being run
#' 
#' 
get_results = function(method,seed,resultsdir="results"){
  data.frame(method=method$name,seed=seed,ldply(seed,get_one_result,method=method,resultsdir=resultsdir))  
}


#' @title Aggregate the results of multiple methods for multiple trials
#'
#' @description Aggregate the results of multiple methods for multiple trials
#' 
#' @param methodslist a list of methods. Each method is itself a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' @param seed (list of integers) the seeds for the pseudo-rng which identifies/indexes trials
#' @return a data frame of results, with one row for each trial/method combination. The details of the columns will depend on the comparison being run
#' 
#' 
#' # aggregate all the results into a data frame
# seed is a vector of seeds to use
# methods is a list of methods to use
aggregate_results = function(methodslist,seed,resultsdir="results"){
  ldply(methodslist,get_results,seed=seed,resultsdir=resultsdir)
}


