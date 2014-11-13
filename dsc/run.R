source("dsc.R")

source("parammaker.R")
source("datamaker.R")
source("score.R")
library("plyr")
library("ashr")


#read in methods
methods = system("ls methods",intern=TRUE)
for(m in 1:length(methods)){
  source(file.path("methods",methods[m]))
}

seedA = 1:10
seedB= 101:111
scenario_seedlist = list(A=seedA,B=seedB)
methods=list()
ashflavorlist = list(hu=list(mixcompdist="halfunif",method="fdr")
                     ,u=list(mixcompdist="uniform",method="fdr"))
methods$ash=list(name="ash",fn="ash.wrapper",flavorlist = ashflavorlist)
methods$qval = list(name="qvalue",fn="qvalue.wrapper")




make_directories(methods,names(scenario_seedlist))
make_params(parammaker,scenario_seedlist)
make_data(datamaker,scenario_seedlist)


apply_method(scenario_seedlist, methods$ash,"hu")
apply_method(scenario_seedlist, methods$ash,"u")
apply_method(scenario_seedlist, methods$qval)

score_method(scenario_seedlist, methods$ash,score,"hu")
score_method(scenario_seedlist, methods$ash,score, "u")
score_method(scenario_seedlist, methods$qval,score)

res=aggregate_results(methods,scenario_seedlist)
aggregate(td~method+flavor+scenario,res,mean)
aggregate(fd~method+flavor+scenario,res,mean)
