#EM algorithm for estimating mixing proportions given TI table,
#parent TI table and beta parameters  
binash=function(tit,ptit,qq,mode=1,maxit=2000,tol=1e-5){
  n=dim(tit)[2]
  J=dim(tit)[1]-1
  nt=tit[-1,]
  ns=ptit[,((1:(2*n))%%2==1)]
  nf=nt-ns
  maxit=maxit
  tol=tol
  convcrit=10.0
  iter=0
  lik=NULL
  likk=0
  if(mode==1){
    pip1=0.5*rep(1,J)
    while((convcrit > tol) & (iter < maxit)){ 
      lik.old=likk
      pi=pip1
      pmat0=log(pi)+nt*log(1/2)
      pmat1=log(1-pi)+lbeta(ns+qq,nf+qq)-lbeta(qq,qq)
      pm=pmax(pmat0,pmat1)
      pmat0=pmat0-pm
      pmat1=pmat1-pm
      pmat=exp(pmat0)/(exp(pmat0)+exp(pmat1))
      #likk=apply(log(exp(pmat0)+exp(pmat1)),1,sum)
      likk=apply(log(exp(pmat0)+exp(pmat1))+pm,1,sum)
      lik=cbind(lik,likk)
      pip1=rowMeans(pmat)
      pip1=pmax(pip1,1e-8)
      p=pip1
      convcrit=max((likk-lik.old)^2)
      iter=iter+1
    }
  }else if(mode==2){
    K=length(qq)
    pip1=1/K*matrix(1,nrow=J,ncol=K)
    qq=qq%o%rep(1,n)
    rK=rep(1,K)
    pmat_i=list(0)
    for(j in 1:J){
      nns=rK%o%ns[j,]
      nnf=rK%o%nf[j,]
      pmat_i[[j]]=lbeta(nns+qq,nnf+qq)-lbeta(qq,qq)
    }
    while((convcrit > tol) & (iter < maxit)){    
      pi=pip1
      lik.old=likk
      for(j in 1:J){
        pp=pi[j,]%o%rep(1,n)
        pmat=log(pp)+pmat_i[[j]]
        pm=apply(pmat,2,max)
        pmat=pmat-rK%o%pm
        pip1.norm=exp(pmat)/(rK%o%colSums(exp(pmat)))
        pip1[j,]=rowMeans(pip1.norm)
        likk[j]=sum(log(colSums(exp(pmat)))+pm)
      }
      lik=cbind(lik,likk)
      pip1=pmax(pip1,1e-8)
      p=pip1
      convcrit=max((likk-lik.old)^2)
      iter=iter+1
    }
  }
  return(list(p=p,lik=lik))
}
