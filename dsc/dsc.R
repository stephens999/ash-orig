library(plyr)
library(reshape2)
#' @title return the path to a data file, parameter file, output file or results file
#'
#' @description return the path to a data file, parameter file, output file or results file
#' 
#' @param indexlist list that ``indexes" the comparison being used, with components
#' \item{seed}{the seed for the pseudo-rng which identifies/indexes the file}
#' \item{scenario}{string indicating scenario name}
#' @param datadir/paramdir/outputdir/resultsdir the (relative) path to the directory containing the relevant files
#' 
#' @return string containing path to file
#' 
#' @example datafile(1)
#' 

datafile = function(indexlist,datadir="data"){
  return(file.path(datadir,indexlist$scenario,paste0("data.",indexlist$seed,".RData")))
}

paramfile = function(indexlist,paramdir="param"){
  return(file.path(paramdir,indexlist$scenario,paste0("param.",indexlist$seed,".RData")))
}

outputfile = function(methodname,indexlist,flavor=NULL, outputdir="output"){
  methodname=long_methodname(methodname,flavor)
  return(file.path(outputdir,indexlist$scenario,methodname,paste0("output.",indexlist$seed,".RData")))
}

resultsfile = function(methodname,indexlist,flavor=NULL, resultsdir="results"){
  methodname=long_methodname(methodname,flavor)    
  return(file.path(resultsdir,indexlist$scenario,methodname,paste0("results.",indexlist$seed,".RData")))
}

#' @title combine a method name and flavor to produce a new method name 
#'
#' @description combine a method name and flavor to produce a new method name 
#' 
#' @param methodname (string) name of a method
#' @param flavor (sting) name of a flavor
#' 
#' @return string method.flavor, or if flavor is null then just method
#' 
long_methodname=function(methodname,flavor=NULL){
  if(is.null(flavor)){ 
    return(methodname)
  }
  else{
    return(paste(methodname,flavor,sep="."))
  }
}

#' @title Apply a method to an input and produce output
#'
#' @description Apply a single method to a single input trial and produce (and save) corresponding output
#' 
#' @param indexlist list that ``indexes" the comparison being used, with components
#' \item{seed}{the seed for the pseudo-rng which identifies/indexes the file}
#' \item{scenario}{string indicating scenario name}
#' 
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method} 
#' @flavor a string indicating which element of methods$flavorlist to use as additional arguments 
#' 
#' @return output a list of appropriate format to be determined by the comparison being run
#' 
#' 
apply_method_singletrial=function(indexlist,method,flavor=NULL){
  load(file=datafile(indexlist))
  if(!is.null(flavor)){
    output=do.call(method$fn,list(input=data$input,add.args=method$flavorlist$flavor))
  } else {
    output=do.call(method$fn,list(input=data$input))
  }
  save(output,file=outputfile(method$name,indexlist,flavor=flavor))
  return(output)
}

#' @title Apply a method to all inputs and produce corresponding outputs
#'
#' @description Apply a method to all inputs and produce corresponding outputs (by repeated application of apply_method_once).
#' Results are saved in a file
#'  
#' @param scenario_seedlist (list), one element for each scenario, each element contains the vector of seeds for the pseudo-rng which identifies/indexes trials
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' 
#' @return none - the results are saved in a file determined by outputfile(seed, method$methodname)
#' 
apply_method = function(scenario_seedlist,method,flavor=NULL){
  combo = melt(scenario_seedlist,value.name="seed")
  names(combo)[2]="scenario"  
  d_ply(combo,.(seed,scenario),apply_method_singletrial,method=method,flavor=flavor)
}

#' @title Score a method on a single trial and save results
#'
#' @description Score results of a single method for a single trial and produce (and save) corresponding results
#' 
#' @param indexlist list that ``indexes" the comparison being used, with components
#' \item{seed}{the seed for the pseudo-rng which identifies/indexes the file}
#' \item{scenario}{string indicating scenario name}
#' 
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' @param scorefn a function that scores output based on comparisons with input, parameters and metadata
#'
#' @return results, a list of appropriate format to be determined by the comparison being run (maybe required to be a dataframe?)
#' 
#' 
score_method_singletrial = function(indexlist,method,scorefn,flavor=NULL){
  load(file=paramfile(indexlist))
  load(file=datafile(indexlist))
  load(file=outputfile(method$name,indexlist,flavor))
  
  results=scorefn(param,data,output)
  save(results,file=resultsfile(method$name,indexlist,flavor))
  return(results)
}

#' @title Score a method on all trials and save results
#'
#' @description Score a method on all trials and save results (by repeated application of score_method_singletrial)
#' 
#' @param scenario_seedlist (list) the seeds for the pseudo-rng for each scenario which identifies/indexes trials
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' @param scorefn a function that scores output based on comparisons with input, parameters and metadata
#'
#' @return results, a list of appropriate format to be determined by the comparison being run (maybe required to be a dataframe?)
#' 
#'
score_method=function(scenario_seedlist,method,scorefn,flavor=NULL){
  combo = melt(scenario_seedlist,value.name="seed")
  names(combo)[2]="scenario"
  d_ply(combo,.(seed,scenario),score_method_singletrial,method=method,scorefn=scorefn,flavor=flavor)
}





#' @title Get the results of a single method for a single trial
#'
#' @description Get the results of a single method for a single trial
#' 
#' @param indexlist list that ``indexes" the comparison being used, with components
#' \item{seed}{the seed for the pseudo-rng which identifies/indexes the file}
#' \item{scenario}{string indicating scenario name}
#' 
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#'
#' @return results, a data frame of results, the details will depend on the comparison being run
#' 
get_one_result = function(indexlist,method,flavor=NULL){
  load(file=resultsfile(method$name,indexlist,flavor))
  return(data.frame(results))
}

#' @title Get the results of a single method for multiple trials (one flavor)
#'
#' @description Get the results of a single method for multiple trials (one flavor)
#' 
#' @param method a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' @param seed (list of integers) the seeds for the pseudo-rng which identifies/indexes trials
#' @return a data frame of results, with one row for each trial. The details of the columns will depend on the comparison being run
#' 
#' 
get_results = function(method,scenario_seedlist,flavor=NULL){
  if(is.null(flavor)){
    flavorname="NA"
  } else {
    flavorname=flavor
  }
  combo = melt(scenario_seedlist,value.name="seed")
  names(combo)[2]="scenario"
    
  data.frame(method=method$name,flavor=flavorname,ddply(combo,.(seed,scenario),get_one_result,method=method,flavor=flavor))  
}

get_results_all_flavors = function(method,scenario_seedlist){
  if(is.null(method$flavorlist)){
    return(get_results(method,scenario_seedlist))
  } else {
    return(ldply(names(method$flavorlist),get_results,method=method,scenario_seedlist=scenario_seedlist))
  }  
}

#' @title Aggregate the results of multiple methods for multiple trials
#'
#' @description Aggregate the results of multiple methods for multiple trials
#' 
#' @param methodslist a list of methods. Each method is itself a list with elements
#' \item{methodname}{string by which method should be identified}
#' \item{methodfn}{name of function that is used to call the method}  
#' \item{flavorlist}{list of flavors of the method (each flavor is a list containing additional parameters to be passed to methodfn)}
#' @param seed (list of integers) the seeds for the pseudo-rng which identifies/indexes trials
#' @return a data frame of results, with one row for each trial/method combination. The details of the columns will depend on the comparison being run
#' 
#' 
#' 
# aggregate all the results into a data frame
# seed is a vector of seeds to use
# methods is a list of methods to use
aggregate_results = function(methodslist,scenario_seedlist){
  ldply(methodslist,get_results_all_flavors,scenario_seedlist=scenario_seedlist)
}


make_directories_singlemethod_singleflavor = function(method,flavor=NULL,scenario=NULL){
  if(is.null(scenario)) scenario = "default_scenario" 
  system(paste0("mkdir ",file.path("output",scenario,long_methodname(method$name,flavor))))
  system(paste0("mkdir ",file.path("results",scenario,long_methodname(method$name,flavor))))
}

 
make_directories_singlemethod = function(method,scenario=NULL){
  if(is.null(method$flavorlist)){
    make_directories_singlemethod_singleflavor(method=method,scenario=scenario)  
  } else {
    lapply(names(method$flavorlist),make_directories_singlemethod_singleflavor,method=method,scenario=scenario)
  }
}

make_scenario_subdirectories_singlescenario=function(methodslist,scenario=NULL){
  if(is.null(scenario)) scenario = "default_scenario" 
  system(paste0("mkdir ",file.path("param",scenario)))
  system(paste0("mkdir ",file.path("data",scenario)))
  system(paste0("mkdir ",file.path("output",scenario)))
  system(paste0("mkdir ",file.path("results",scenario)))
  l_ply(methodslist,make_directories_singlemethod,scenario)
}

make_directories = function(methodslist,scenario=NULL){
  system("mkdir param")
  system("mkdir data")
  system("mkdir output")
  system("mkdir results")
  if(is.null(scenario)) scenario="default_scenario"
  l_ply(scenario,make_scenario_subdirectories_singlescenario,methodslist=methodslist)
}


#' @title Make all the parameters for a DSC
#'
#' @description Make all the parameters for DSC using all combinations of seed and scenario
#' 
#' @param parammaker a function for making parameters from seeds and scenario combinations
#' @param scenario_seedlist a list of integer vectors, one for each scenario. Each list gives the seeds to be used for that scenario.
#' 
#' @return parameters are saved in files in the params subdirectory
#' 
make_params = function(parammaker,scenario_seedlist){
  combo = melt(scenario_seedlist,value.name="seed")
  names(combo)[2]="scenario"
  d_ply(combo,.(seed,scenario),parammaker)
}


#' @title Make all the inputs for a DSC
#'
#' @description Make all the parameters for DSC using all combinations of seed and scenario
#' 
#' @param parammaker a function for making parameters from seeds and scenario combinations
#' @param seed (list of integers) the seeds for the pseudo-rng which identifies/indexes trials
#' @param scenario (list of strings) the names of the scenarios 
#' 
#' @return parameters are saved in files in the params subdirectory
#' 
make_data = function(datamaker,scenario_seedlist){
  combo = melt(scenario_seedlist,value.name="seed")
  names(combo)[2]="scenario"
  d_ply(combo,.(seed,scenario),datamaker)
}


