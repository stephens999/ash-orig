library(dscr)
source_dir("methods")

###### Initialize #######

dsc_robust=new_dsc("robust","dsc-robust-files")

###### Add Scenarios #####
source_dir("../dsc-shrink/datamakers")
source("../dsc-shrink/add_named_scenarios.R")
add_named_scenarios(dsc_robust,c("spiky","near-normal","flat-top","skew","big-normal","bimodal"))


###### Add Methods #####

add_method(dsc_robust,"ash.u",ash.wrapper,args=list(mixcompdist="uniform",cxx=TRUE))
add_method(dsc_robust,"ash.hu",ash.wrapper,args=list(mixcompdist="halfuniform",cxx=TRUE))
add_method(dsc_robust,"ash.n",ash.wrapper,args=list(mixcompdist="normal",cxx=TRUE))



####### Define Score and Add it #######

score = function(data, output){
  x=output$loglik
  #names(x)=paste0('C',1:length(x))
  #class(x)<-'data.frame'
  return(list(diff1 = max(x)-x[1],diff2=(max(x)-min(x))))
}

add_score(dsc_robust,score)

######## Run the DSC #################

res_robust=run_dsc(dsc_robust)
save(dsc_robust,file="dsc-robust-files/dsc_robust.RData")






