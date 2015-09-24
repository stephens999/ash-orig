source_dir("datamakers")


#Now, for each scenario create an element of scenarios of the following form
add_scenario(dsc_shrink,name="A",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(c(2/3,1/3),c(0,0),c(1,2)),
                      min_pi0=0,
                      max_pi0=0,
                      nsamp=100,
                      betahatsd=1
                      ),
                    seed=1:100)

add_scenario(dsc_shrink,name="B",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
                      min_pi0=0,
                      max_pi0=0,
                      nsamp=100,
                      betahatsd=1
                    ),
                    seed=1:100)


add_scenario(dsc_shrink,name="C",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1/4,1/4,1/3,1/6),c(-2,-1,0,1),c(2,1.5,1,1)),
                        min_pi0=0,
                        max_pi0=0,
                        nsamp=100,
                        betahatsd=1
                      ),
                      seed=1:100)

#scenarios An, Bn, Cn are the same as A,B,C but with nulls included (pi0 uniform on [0,1])
add_scenario(dsc_shrink,name="An",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(c(2/3,1/3),c(0,0),c(1,2)),
                      min_pi0=0,
                      max_pi0=1,
                      nsamp=100,
                      betahatsd=1
                    ),
                    seed=1:100)

add_scenario(dsc_shrink,name="Bn",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
                      min_pi0=0,
                      max_pi0=1,
                      nsamp=100,
                      betahatsd=1
                    ),
                    seed=1:100)


add_scenario(dsc_shrink,name="Cn",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1/4,1/4,1/3,1/6),c(-2,-1,0,1),c(2,1.5,1,1)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=100,
                        betahatsd=1
                      ),
                      seed=1:100)

add_scenario(dsc_shrink,name="hard",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(.4,.2,.2,.2),c(0,0,0,0),c(.25,.5,1,2)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=100,
                        betahatsd=1
                      ),
                      seed=1:100)

add_scenario(dsc_shrink,name="easy",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1),c(0),c(4)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=100,
                        betahatsd=1
                      ),
                      seed=1:100)

add_scenario(dsc_shrink,name="efron_FCR",
            fn=rnormmix_datamaker,
            args=list(
              g=normalmix(c(1),c(-3),c(1)),
              min_pi0=0.9,
              max_pi0=0.9,
              nsamp=100,
              betahatsd=1
            ),
            seed=1:100)



