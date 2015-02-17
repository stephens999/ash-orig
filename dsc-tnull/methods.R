sourceDir("methods")
methods=list()
methods[[1]] = list(name="ash.hu",fn =ash.wrapper,args=list(mixcompdist="halfunif",method="fdr"))
methods[[2]] = list(name="ash.hu.betaeff",fn =ash.betaeff.wrapper,args=list(mixcompdist="halfunif",method="fdr"))
methods[[3]] = list(name="ash.hu.sigmaeff",fn =ash.sigmaeff.wrapper,args=list(mixcompdist="halfunif",method="fdr"))
methods[[4]] = list(name="qvalue",fn =qvalue.wrapper)
methods[[5]] = list(name="ash.hu.ES",fn =ash.wrapper,args=list(mixcompdist="halfunif",method="fdr",model="ES"))


#methods[[2]] = list(name="ash.u",fn =ash.wrapper,args=list(mixcompdist="unif",method="fdr"))
#methods[[3]] = list(name="ash.n",fn =ash.wrapper,args=list(mixcompdist="normal",method="fdr"))
