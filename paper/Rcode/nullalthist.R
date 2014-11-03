#plot a histogram of z scores, highlighting the alternative distribution
#of z scores that is implied by localfdr values lfdr.
nullalthist = function(z,lfdr,nullcol="blue",altcol="cyan",ncz=100,...){
  h=hist(z, freq=FALSE,col=nullcol,nclass=ncz,...)
  avlfdr = unlist(lapply(split(lfdr,cut(z,h$breaks),drop=FALSE),mean))
  h$density = (1-avlfdr) * h$density
  plot(h,add=TRUE,col=altcol,freq=FALSE)
}

#this one puts the null on the bottom
altnullhist = function(z,lfdr,nullcol="blue",altcol="cyan",ncz=100,...){
  h=hist(z, freq=FALSE,col=altcol,nclass=ncz,...)
  avlfdr = unlist(lapply(split(lfdr,cut(z,h$breaks),drop=FALSE),mean))
  h$density = avlfdr * h$density
  plot(h,add=TRUE,col=nullcol,freq=FALSE)
}

plotall_hist=function(sim,iter=1,histfun=nullalthist){
  hh.zscore=sim$zscore[[iter]]    
  par(mfcol=c(2,2))
  histfun(hh.zscore,sim$betahat.fdrtool[[iter]]$lfdr,main="fdrtool")  
  histfun(hh.zscore,sim$betahat.locfdr[[iter]]$fdr,main="locfdr")
  histfun(hh.zscore,sim$betahat.mixfdr[[iter]]$fdr,main="mixfdr")
  histfun(hh.zscore,sim$betahat.ash.n[[iter]]$lfdr,main="ash")
  par(mfcol=c(1,1))
}
