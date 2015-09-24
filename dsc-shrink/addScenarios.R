source_dir("datamakers")


addNamedScenarios = function(dsc,names,nsamp=1000){

#the -nn suffix indicates no nulls (min_pi0 = max_pi0 = 0)

if("near-normal-nn" %in% names){
	add_scenario(dsc,name="near-normal-nn",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(c(2/3,1/3),c(0,0),c(1,2)),
                      min_pi0=0,
                      max_pi0=0,
                      nsamp=nsamp,
                      betahatsd=1
                      ),
                    seed=1:100)
}

if("flat-top-nn" %in% names){
	add_scenario(dsc,name="flat-top-nn",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
                      min_pi0=0,
                      max_pi0=0,
                      nsamp=nsamp,
                      betahatsd=1
                    ),
                    seed=1:100)
}

add_scenario(dsc,name="skew-nn",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1/4,1/4,1/3,1/6),c(-2,-1,0,1),c(2,1.5,1,1)),
                        min_pi0=0,
                        max_pi0=0,
                        nsamp=nsamp,
                        betahatsd=1
                      ),
                      seed=1:100)

add_scenario(dsc,name="near-normal",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(c(2/3,1/3),c(0,0),c(1,2)),
                      min_pi0=0,
                      max_pi0=1,
                      nsamp=nsamp,
                      betahatsd=1
                    ),
                    seed=1:100)

add_scenario(dsc,name="flat-top",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
                      min_pi0=0,
                      max_pi0=1,
                      nsamp=nsamp,
                      betahatsd=1
                    ),
                    seed=1:100)


add_scenario(dsc,name="skew",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1/4,1/4,1/3,1/6),c(-2,-1,0,1),c(2,1.5,1,1)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=nsamp,
                        betahatsd=1
                      ),
                      seed=1:100)

add_scenario(dsc,name="spiky",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(.4,.2,.2,.2),c(0,0,0,0),c(.25,.5,1,2)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=nsamp,
                        betahatsd=1
                      ),
                      seed=1:100)

add_scenario(dsc,name="big-normal",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1),c(0),c(4)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=nsamp,
                        betahatsd=1
                      ),
                      seed=1:100)


add_scenario(dsc,name="bimodal",
            fn=rnormmix_datamaker,
            args=list(
              g=normalmix(c(0.5,0.5),c(-2,2),c(1,1)),
              min_pi0=0,
              max_pi0=1,
              nsamp=nsamp,
              betahatsd=1
            ),
            seed=1:100)


}
