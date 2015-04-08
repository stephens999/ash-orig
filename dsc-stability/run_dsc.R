library(dscr)
sourceDir("methods")
sourceDir("datamakers")

###### Initialize #######

dsc_ar=new.dsc("ash-robust","ash-robust-files")

###### Add Scenarios #####


#Now, for each scenario create an element of scenarios of the following form
addScenario(dsc_ar,name="A",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(c(2/3,1/3),c(0,0),c(1,2)),
                      min_pi0=0,
                      max_pi0=0,
                      nsamp=1000,
                      betahatsd=1
                    ),
                    seed=1:100)

addScenario(dsc_ar,name="B",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
                      min_pi0=0,
                      max_pi0=0,
                      nsamp=1000,
                      betahatsd=1
                    ),
                    seed=1:100)


addScenario(dsc_ar,name="C",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1/4,1/4,1/3,1/6),c(-2,-1,0,1),c(2,1.5,1,1)),
                        min_pi0=0,
                        max_pi0=0,
                        nsamp=1000,
                        betahatsd=1
                      ),
                      seed=1:100)

#scenarios An, Bn, Cn are the same as A,B,C but with nulls included (pi0 uniform on [0,1])
addScenario(dsc_ar,name="An",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(c(2/3,1/3),c(0,0),c(1,2)),
                      min_pi0=0,
                      max_pi0=1,
                      nsamp=1000,
                      betahatsd=1
                    ),
                    seed=1:100)

addScenario(dsc_ar,name="Bn",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
                      min_pi0=0,
                      max_pi0=1,
                      nsamp=1000,
                      betahatsd=1
                    ),
                    seed=1:100)


addScenario(dsc_ar,name="Cn",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1/4,1/4,1/3,1/6),c(-2,-1,0,1),c(2,1.5,1,1)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=1000,
                        betahatsd=1
                      ),
                      seed=1:100)

addScenario(dsc_ar,name="hard",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(.4,.2,.2,.2),c(0,0,0,0),c(.25,.5,1,2)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=1000,
                        betahatsd=1
                      ),
                      seed=1:100)

addScenario(dsc_ar,name="easy",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1),c(0),c(4)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=1000,
                        betahatsd=1
                      ),
                      seed=1:100)


addScenario(dsc_ar,name="bimodal",
            fn=rnormmix_datamaker,
            args=list(
              g=normalmix(c(0.5,0.5),c(-2,2),c(1,1)),
              min_pi0=0,
              max_pi0=1,
              nsamp=1000,
              betahatsd=1
            ),
            seed=1:100)



# addScenario(dsc_ar,name="hard-b",
#                       fn=rnormmix_datamaker,
#                       args=list(
#                         g=normalmix(c(.4,.2,.2,.2),c(0,0,0,0),c(.25,.5,1,2)),
#                         min_pi0=0,
#                         max_pi0=1,
#                         nsamp=10000,
#                         betahatsd=1
#                       ),
#                       seed=1:100)





###### Add Methods #####

addMethod(dsc_ar,"ash.u",ash.wrapper,args=list(mixcompdist="uniform",cxx=TRUE))
addMethod(dsc_ar,"ash.hu",ash.wrapper,args=list(mixcompdist="halfuniform",cxx=TRUE))
addMethod(dsc_ar,"ash.n",ash.wrapper,args=list(mixcompdist="normal",cxx=TRUE))



####### Define Score and Add it #######

score = function(data, output){
  x=output$loglik
  #names(x)=paste0('C',1:length(x))
  #class(x)<-'data.frame'
  return(list(diff1 = max(x)-x[1],diff2=(max(x)-min(x))))
}

addScore(dsc_ar,score)

######## Run the DSC #################

#reset_dsc(dsc_ar,force=TRUE)
res=run_dsc(dsc_ar)
save(dsc_ar,file="dsc_ar.RData")

quantile((res %>% filter(method=="ash.n"))$diff1,c(0.95,0.99,1))
quantile((res %>% filter(method=="ash.hu"))$diff1,c(0.8,0.9,0.95,0.99,1))
quantile((res %>% filter(method=="ash.u"))$diff1,c(0.8,0.9,0.95,0.99,1))
plot(ecdf((res %>% filter(method %in% c("ash.u","ash.hu")))$diff1))
quantile((res %>% filter(method %in% c("ash.u","ash.hu")))$diff1,c(0.8,0.9,0.95,0.99,1))
mean((res %>% filter(method %in% c("ash.u","ash.hu")))$diff1>1)



plot(ecdf((res %>% filter(method %in% c("ash.n")))$diff1))
        




