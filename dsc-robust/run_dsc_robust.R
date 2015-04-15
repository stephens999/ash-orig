library(dscr)
sourceDir("methods")

###### Initialize #######

dsc_robust=new.dsc("robust","dsc-robust-files")

###### Add Scenarios #####
sourceDir("../dsc-shrink/datamakers")
source("../dsc-shrink/addNamedScenarios.R")
addNamedScenarios(dsc_robust,c("spiky","near-normal","flat-top","skew","big-normal","bimodal"))


###### Add Methods #####

addMethod(dsc_robust,"ash.u",ash.wrapper,args=list(mixcompdist="uniform",cxx=TRUE))
addMethod(dsc_robust,"ash.hu",ash.wrapper,args=list(mixcompdist="halfuniform",cxx=TRUE))
addMethod(dsc_robust,"ash.n",ash.wrapper,args=list(mixcompdist="normal",cxx=TRUE))



####### Define Score and Add it #######

score = function(data, output){
  x=output$loglik
  #names(x)=paste0('C',1:length(x))
  #class(x)<-'data.frame'
  return(list(diff1 = max(x)-x[1],diff2=(max(x)-min(x))))
}

addScore(dsc_robust,score)

######## Run the DSC #################

#reset_dsc(dsc_robust,force=TRUE)
res_robust=run_dsc(dsc_robust)
#save(res_robust,file="res_robust.RData")
save(dsc_robust,file="dsc_robust.RData")






