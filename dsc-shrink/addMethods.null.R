source_dir("methods")

add_method(dsc_znull,name="ash.hu",fn =ash.wrapper,args=list(mixcompdist="halfunif"),outputtype = "ash_output")
add_method(dsc_znull,name="ash.u",fn =ash.wrapper,args=list(mixcompdist="unif"),outputtype = "ash_output")
add_method(dsc_znull,name="ash.n",fn =ash.wrapper,args=list(mixcompdist="normal"),outputtype = "ash_output")

add_method(dsc_znull,name="ash.hu.s",fn =ash.wrapper,args=list(mixcompdist="halfunif",method="shrink"),outputtype = "ash_output")
add_method(dsc_znull,name="ash.u.s",fn =ash.wrapper,args=list(mixcompdist="unif",method="shrink"),outputtype = "ash_output")
add_method(dsc_znull,name="ash.n.s",fn =ash.wrapper,args=list(mixcompdist="normal",method="shrink"),outputtype = "ash_output")

add_method(dsc_znull,name="ash.hu.nw2",fn =ash.wrapper,args=list(mixcompdist="halfunif",nullweight=2),outputtype = "ash_output")
add_method(dsc_znull,name="ash.u.nw2",fn =ash.wrapper,args=list(mixcompdist="unif",nullweight=2),outputtype = "ash_output")
add_method(dsc_znull,name="ash.n.nw2",fn =ash.wrapper,args=list(mixcompdist="normal",nullweight=2),outputtype = "ash_output")


add_method(dsc_znull,name="ash.hu.s.gmfine",fn =ash.wrapper,args=list(mixcompdist="halfunif",method="shrink",gridmult=2^.25),outputtype = "ash_output")
add_method(dsc_znull,name="ash.u.s.gmfine",fn =ash.wrapper,args=list(mixcompdist="unif",method="shrink",gridmult=2^.25),outputtype = "ash_output")
add_method(dsc_znull,name="ash.n.s.gmfine",fn =ash.wrapper,args=list(mixcompdist="normal",method="shrink",gridmult=2^.25),outputtype = "ash_output")


