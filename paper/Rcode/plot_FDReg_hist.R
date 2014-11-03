plot_FDReg_hist = function(hh.pval,pi0,nc=40,nullcol="blue",altcol="cyan",type=4,textsize=1.2,title="Distribution of p values",...){
  hh.hist=hist(hh.pval,freq=FALSE,xlab="p value",main=title,nclass=nc,col=altcol,...)
  if(type>1){
    abline(h=pi0,col=nullcol,lwd=2)
    
    hh.hist$density=rep(pi0,length(hh.hist$density))  
    #hh.hist$counts=rep(hh.q$pi0*length(hh.pval)/nc,length(hh.hist$counts)) 
    plot(hh.hist,add=TRUE,col=nullcol,freq=FALSE)
  }
  if(type>2){
    abline(v=0.1,lwd=2,col=2)
  }
  if(type>3){
    text(0.05,1.2,labels="A",col=2,cex=1.2)  
    text(0.05,0.4,labels="B",col=2,cex=1.2)  
    text(0.6,3,labels=paste0("FDR = B/(A+B) =  ",round(pi0*0.1*length(hh.pval)/sum(hh.pval<0.1),2)),cex=textsize)
  }
}


