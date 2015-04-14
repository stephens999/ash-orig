sourceDir("methods")

addMethod(dsc_shrink,name="ash.hu.nocxx",fn =ash.wrapper,args=list(mixcompdist="halfunif",cxx=FALSE),outputtype = "ash_output")
addMethod(dsc_shrink,name="ash.hu",fn =ash.wrapper,args=list(mixcompdist="halfunif",cxx=TRUE),outputtype = "ash_output")
addMethod(dsc_shrink,name="ash.u",fn =ash.wrapper,args=list(mixcompdist="unif",cxx=TRUE),outputtype = "ash_output")
addMethod(dsc_shrink,name="ash.n",fn =ash.wrapper,args=list(mixcompdist="normal",cxx=TRUE),outputtype = "ash_output")

addMethod(dsc_shrink,name="ash.hu.nw2",fn =ash.wrapper,args=list(mixcompdist="halfunif",nullweight=2,cxx=TRUE),outputtype = "ash_output")
addMethod(dsc_shrink,name="ash.u.nw2",fn =ash.wrapper,args=list(mixcompdist="unif",nullweight=2,cxx=TRUE),outputtype = "ash_output")
addMethod(dsc_shrink,name="ash.n.nw2",fn =ash.wrapper,args=list(mixcompdist="normal",nullweight=2,cxx=TRUE),outputtype = "ash_output")


addMethod(dsc_shrink,name="ash.hu.s",fn =ash.wrapper,args=list(mixcompdist="halfunif",method="shrink",cxx=TRUE),outputtype = "ash_output")
addMethod(dsc_shrink,name="ash.u.s",fn =ash.wrapper,args=list(mixcompdist="unif",method="shrink",cxx=TRUE),outputtype = "ash_output")
addMethod(dsc_shrink,name="ash.n.s",fn =ash.wrapper,args=list(mixcompdist="normal",method="shrink",cxx=TRUE),outputtype = "ash_output")


addMethod(dsc_shrink,name="mixfdr.tnull", fn=mixfdr.wrapper, args = list(theonull=TRUE),outputtype = "mixfdr_output")
addMethod(dsc_shrink,name="mixfdr.enull", fn=mixfdr.wrapper, args = list(theonull=FALSE),outputtype = "mixfdr_output")

addMethod(dsc_shrink,name="locfdr", fn=locfdr.wrapper,outputtype = "locfdr_output")
addMethod(dsc_shrink,name="qvalue", fn=qvalue.wrapper,outputtype = "qvalue_output")

#addMethod(dsc_shrink,name="mixfdr.tnull.J10", fn=mixfdr.wrapper, args = list(theonull=TRUE,J=10),outputtype = "mixfdr_output")
#addMethod(dsc_shrink,name="mixfdr.enull.J10", fn=mixfdr.wrapper, args = list(theonull=FALSE,J=10),outputtype = "mixfdr_output")
#addMethod(dsc_shrink,name="mixfdr.tnull.J10P0", fn=mixfdr.wrapper, args = list(theonull=TRUE,J=10,P=0),outputtype = "mixfdr_output")
