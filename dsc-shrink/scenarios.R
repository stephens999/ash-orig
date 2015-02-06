sourceDir("datamakers")
scenarios=list()

#Now, for each scenario create an element of scenarios of the following form
scenarios[[1]]=list(name="A",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(c(2/3,1/3),c(0,0),c(1,2)),
                      min_pi0=0,
                      max_pi0=0,
                      nsamp=1000,
                      betahatsd=1
                      ),
                    seed=1:100)

scenarios[[2]]=list(name="B",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
                      min_pi0=0,
                      max_pi0=0,
                      nsamp=1000,
                      betahatsd=1
                    ),
                    seed=1:100)


scenarios[[3]] = list(name="C",
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
scenarios[[4]]=list(name="An",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(c(2/3,1/3),c(0,0),c(1,2)),
                      min_pi0=0,
                      max_pi0=1,
                      nsamp=1000,
                      betahatsd=1
                    ),
                    seed=1:100)

scenarios[[5]]=list(name="Bn",
                    fn=rnormmix_datamaker,
                    args=list(
                      g=normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7)),
                      min_pi0=0,
                      max_pi0=1,
                      nsamp=1000,
                      betahatsd=1
                    ),
                    seed=1:100)


scenarios[[6]] = list(name="Cn",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1/4,1/4,1/3,1/6),c(-2,-1,0,1),c(2,1.5,1,1)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=1000,
                        betahatsd=1
                      ),
                      seed=1:100)

scenarios[[7]] = list(name="hard",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(.4,.2,.2,.2),c(0,0,0,0),c(.25,.5,1,2)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=1000,
                        betahatsd=1
                      ),
                      seed=1:100)

scenarios[[8]] = list(name="easy",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(1),c(0),c(4)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=1000,
                        betahatsd=1
                      ),
                      seed=1:100)



scenarios[[9]] = list(name="hard-b",
                      fn=rnormmix_datamaker,
                      args=list(
                        g=normalmix(c(.4,.2,.2,.2),c(0,0,0,0),c(.25,.5,1,2)),
                        min_pi0=0,
                        max_pi0=1,
                        nsamp=10000,
                        betahatsd=1
                      ),
                      seed=1:100)


#scenarios[[2]]=list(name="2",fn=datamaker,args=list(g=normalmix(c(0.98,0.02),c(0,0),c(0,3)),nsamp=1000,betahatsd=1),seed=1:100)
#scenarios[[3]]=list(name="3",fn=datamaker,args=list(g=normalmix(c(0.98,0.02),c(0,0),c(0,3)),nsamp=1000,betahatsd=1),seed=1:100)
#scenarios[[4]]]=list(name="1b",fn=datamaker,args=list(g=normalmix(c(0.5,0.5),c(0,0),c(0,3)),nsamp=10000,betahatsd=1),seed=1:100)

