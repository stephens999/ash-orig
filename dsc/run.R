source("dsc.R")

source("parammaker.R")
source("datamaker.R")
source("score.R")
library("plyr")

system("mkdir param")
system("mkdir data")
system("mkdir output")
system("mkdir output/ash")
system("mkdir output/qvalue")
system("mkdir results")
system("mkdir results/ash")
system("mkdir results/qvalue")

#read in methods
methods = system("ls methods",intern=TRUE)
for(m in 1:length(methods)){
  source(file.path("methods",methods[m]))
}

seed = 1:50

l_ply(seed,parammaker)
l_ply(seed,datamaker)

ashmethod=list(name="ash",fn="ash.wrapper")
qvalmethod = list(name="qvalue",fn="qvalue.wrapper")

apply_method(seed, ashmethod)
apply_method(seed, qvalmethod)

score_method(seed, ashmethod,score)
score_method(seed, qvalmethod,score)

aggregate_results(list(ashmethod,qvalmethod),seed)
