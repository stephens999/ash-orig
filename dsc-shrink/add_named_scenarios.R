source_dir("datamakers")

gdef = list(spiky=normalmix(c(.4,.2,.2,.2),c(0,0,0,0),c(.25,.5,1,2)),
            skew=normalmix(c(1/4,1/4,1/3,1/6),c(-2,-1,0,1),c(2,1.5,1,1)),
            bignormal=normalmix(c(1),c(0),c(4)),
            bimodal=normalmix(c(0.5,0.5),c(-2,2),c(1,1)),
            flat_top=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
            near_normal=normalmix(c(2/3,1/3),c(0,0),c(1,2))
            )

add_named_scenarios = function(dsc,names,nsamp=1000,min_pi0=0,max_pi0=1,suffix=""){

if("near-normal" %in% names){
  add_scenario(dsc,name=paste0("near-normal",suffix),
                    fn=rnormmix_datamaker,
                    args=list(
                      g=gdef$near_normal,
                      min_pi0=min_pi0,
                      max_pi0=max_pi0,
                      nsamp=nsamp,
                      betahatsd=1
                    ),
                    seed=1:100)
}

if("flat-top" %in% names){
  add_scenario(dsc,name=paste0("flat-top",suffix),
                    fn=rnormmix_datamaker,
                    args=list(
                      g=gdef$flat_top,
                      min_pi0=min_pi0,
                      max_pi0=max_pi0,
                      nsamp=nsamp,
                      betahatsd=1
                    ),
                    seed=1:100)
}

if("skew" %in% names){
  add_scenario(dsc,name=paste0("skew",suffix),
                      fn=rnormmix_datamaker,
                      args=list(
                        g=gdef$skew,
                        min_pi0=min_pi0,
                        max_pi0=max_pi0,
                        nsamp=nsamp,
                        betahatsd=1
                      ),
                      seed=1:100)
}

if("spiky" %in% names){
  add_scenario(dsc,name=paste0("spiky",suffix),
                      fn=rnormmix_datamaker,
                      args=list(
                        g=gdef$spiky,
                        min_pi0=min_pi0,
                        max_pi0=max_pi0,
                        nsamp=nsamp,
                        betahatsd=1
                      ),
                      seed=1:100)
}

if("big-normal" %in% names){
  add_scenario(dsc,name=paste0("big-normal",suffix),
                      fn=rnormmix_datamaker,
                      args=list(
                        g=gdef$bignormal,
                        min_pi0=min_pi0,
                        max_pi0=max_pi0,
                        nsamp=nsamp,
                        betahatsd=1
                      ),
                      seed=1:100)
}

if("bimodal" %in% names){
  add_scenario(dsc,name=paste0("bimodal",suffix),
            fn=rnormmix_datamaker,
            args=list(
              g=gdef$bimodal,
              min_pi0=min_pi0,
              max_pi0=max_pi0,
              nsamp=nsamp,
              betahatsd=1
            ),
            seed=1:100)
}

}
