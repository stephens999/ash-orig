scenarios=list()
scenarios[[1]]=list(name="A",fn=datamaker,args=list(g=normalmix(c(0.5,0.5),c(0,0),c(0,3)),nsamp=1000),seed=1:100)
scenarios[[2]]=list(name="B",fn=datamaker,args=list(g=normalmix(c(0.98,0.02),c(0,0),c(0,3)),nsamp=1000),seed=1:100)
