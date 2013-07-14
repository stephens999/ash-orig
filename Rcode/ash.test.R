source("ash.R")

#Test VBEM
abf = rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,1,0),c(0,0,1,0))
eps = 1e-10
abf[abf==0] = eps #replace 0 with small number
print(all.equal(VBEM(abf,c(1,1,1,1))$post,c(2,2,4,1)))
print(all.equal(VBEM(abf,c(1,2,1,1))$post,c(2,3,4,1)))
