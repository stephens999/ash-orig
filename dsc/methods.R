sourceDir("methods")
methods=list()
methods[[1]] = list(name="ash.hu",fn =ash.wrapper,args=list(mixcompdist="halfunif",method="fdr"))
methods[[2]] = list(name="ash.u",fn =ash.wrapper,args=list(mixcompdist="unif",method="fdr"))
methods[[3]] = list(name="ash.n",fn =ash.wrapper,args=list(mixcompdist="normal",method="fdr"))
methods[[4]] = list(name="qvalue",fn=qvalue.wrapper,args=NULL)

