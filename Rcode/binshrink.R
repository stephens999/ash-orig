#This is a R port of the original BMSMShrink function by Kolaczyk, with some modifications and improvements


ash.repodir = scan(".ash.repodir.txt",what=character())
source(file.path(ash.repodir,"/Rcode/bash.R")) 
require(Rcpp)
require(inline)

#interleave two vectors
interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}


rshift = function(x){L=length(x); return(c(x[L],x[-L]))}

lshift = function(x){return(c(x[-1],x[1]))}


#The following produces both a TItable,
#and a "parent" table whose pairwise comparisons would be used to create a TI table
#So, for example, in the ith row, elements 1, 2 would be the parents of the
#first element in the (i+1)the row of the TI table
# INPUT: sig, an n vector of Poisson counts at n locations
# OUTPUT: a list, with elements
# TItable - the usual TI table
# parent - the parent values used to make the TI table
ParentTItable=function(sig){
  n = length(sig)
  J = log2(n)
  
# Create decomposition table of signal, using pairwise sums,
# keeping just the values that are *not* redundant under the
# shift-invariant scheme.  This is very similar to TI-tables
# in Donoho and Coifman's TI-denoising framework.
  dmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  #dmat[1,] = as.matrix(sig)  
  dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table
  
  for(D in 0:(J-1)){
    nD = 2^(J-D); 
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      dmat2[D+1,ind2] = c(x,rx)
    }
  }
  return(list(TItable=dmat,parent=dmat2))
}



#inc <- '#include <cmath>
#        NumericVector rshift(NumericVector x){
#          L=x.size; 
#          double xx=x(L);
#          x.erase(L);
#          x.insert(x.begin(),xx);
#          return(x)
#        }
#        NumericVector lshift(NumericVector x){
#          double xx=x(1);
#          x.erase(1);
#          x.insert(x.end(),xx)
#          return(x)
#        }
#        '


#cxxParentTItable
#cxxParentTItable implements ParentTItable in C++
#note that input to cxxParentTItable is an nsig by n matrix sig
#while input to ParentTItable is a 1 by n vector instead
#src is a string containing the C++ code
src <- '
        NumericVector signal=sig; 
        int n=(int) signal.size();
        int J=(int) log2((double)n);

        NumericMatrix parent(J,2*n);
        NumericMatrix TItable(J+1,n);
        TItable(0,_) = signal;
        for (int D=0; D<J; D++){
           int nD=(int) pow(2., (double) (J-D)), pD=(int) pow(2.,(double) D);
           for (int l=0; l<pD; l++){
              int a=l*nD+1, b=2*l*nD+1;
              double d;
              for (int i=0; i<nD-1; i++){
                 d=TItable(D,a+i-1);
                 parent(D,b+i-1)=d;
                 parent(D,b+i+nD)=d;
              }
              //i=nD-1
              d=TItable(D,a+nD-2);
              parent(D,b+nD-2)=d;
              parent(D,b+nD-1)=d;

              for (int i=0; i<nD; i++)
                TItable(D+1,a+i-1)=parent(D,b+2*i-1)+parent(D,b+2*i);
          }
        }
        return(List::create(Named("TItable")=TItable, Named("parent")=parent));
        '
cxxParentTItable <- cxxfunction(signature(sig="numeric"),
                                body=src,
                                plugin="Rcpp",
                                inc="#include <cmath>")


#This function computes the posterior means
sfunc=function(p,q,x0,x1,mode){
  if(mode==1){
    xx = x0 + x1
    if(p==1){
      ss = 0.5*rep(1,length(x0))
    }else if(p == 0){
      ss = (x0 + q)/(xx + 2*q)
    }else{
# Compute first half denomenator of sum.
      fhd = log(1-p) - log(p) + (xx+1)*log(2) + lbeta(x0+q,x1+q) - lbeta(q,q)
      fhd = 2 + exp(fhd)
     
# Compute second half denomenator.
      shd1 = log(p) - log(1-p) + lbeta(q,q) - xx*log(2) - lbeta(x0+q+1,x1+q)
      shd  = exp(shd1) + (xx + 2*q)/(x0 + q)
     
# Put together two pieces.  Numerators are just 1.
      ss = (1/fhd) + (1/shd)
    }
  }else if(mode==2){
    nq = length(q)
    nn = length(x0)
    x0m = rep(1,nq)%o%x0
    x1m = rep(1,nq)%o%x1
    pm = p%o%rep(1,nn)
    qm = q%o%rep(1,nn)
    rq = rep(1,nq)
    num = log(pm)+lbeta(x0m+qm+1,x1m+qm)-lbeta(qm,qm)
    numpm = apply(num,2,max)
    num = num-rq%o%numpm
    num = log(colSums(exp(num)))+numpm
    den = log(pm)+lbeta(x0m+qm,x1m+qm)-lbeta(qm,qm)
    denpm = apply(den,2,max)
    den = den-rq%o%denpm
    den = log(colSums(exp(den)))+denpm
    ss = exp(num-den)
  }
  return(ss)
}

#This is similar to reverse.pwave used in cyclespin.smooth
reverse.pp=function(dmat,pp,qq,mode){
  n=dim(dmat)[2]
  J=log2(n)  
  if(mode==1){
    qq=rep(qq,J)
  }else if(mode==2){
    qq=rep(1,J)%o%qq
  }
# Beginning with the total number of counts, working from coarse
# scales (i.e., deep depths) upwards, gradually build up the MSPB
# shrinkage estimate by multiplying by appropriate shrinkage factors
# at each level.
  est = dmat[J+1,]
  for (D in J:1){
    nD = 2^(J-D+1)
    nDo2 = nD/2
    for (l in 0:(2^(D-1)-1)){
# Set indexing so as to pick off blocks of size 2^(J-D+1)
# when shrinking estimates at depth D+1 down to finer
# scale at depth D.
      ind = (l*nD+1):((l+1)*nD)
      xxD = dmat[D,ind]
      estvec = est[ind]
    
# In the first half of the vector of D+1-depth estimates,
# we can shrink using the D-depth counts in the order
# in which they appear.
      estl = estvec[1:nDo2]
      xxDl = xxD[seq(1,nD-1,2)]
      xxDr = xxD[seq(2,nD,2)]
      if(mode==1){
        ss = sfunc(pp[D],qq[D],xxDl,xxDr,1)
      }else if(mode==2){
        ss = sfunc(pp[D,],qq[D,],xxDl,xxDr,2)
      }
      nestl = interleave(estl*ss,estl*(1 - ss))
    
# In the second half of the vector of D+1-depth counts,
# we right-shift the D-depth counts, compute the shrunken
# values, and left-shift these values back to the order
# of the above.
      estr = estvec[(nDo2+1):nD]
      sxxD = rshift(xxD)
      xxDl = sxxD[seq(1,nD-1,2)]
      xxDr = sxxD[seq(2,nD,2)]
      if(mode==1){
        ss = sfunc(pp[D],qq[D],xxDl,xxDr,1)
      }else if(mode==2){
        ss = sfunc(pp[D,],qq[D,],xxDl,xxDr,2)
      }
      nestr = interleave(estr*ss,estr*(1 - ss))
      nestr = lshift(nestr)
    
# Combine the estimates from both halves of the D+1-depth
# counts, and store.
      est[ind] = 0.5*( nestl + nestr )
    } 
  } 
  return(est)
}

#The main shrinkage function, as with the original BMSMShrink
#x is the signal
#q is the choice of beta parameters: for mixture of point mass and a single beta distn
#q is the parameter for that beta distn, and for multiple betas q should be a vector of parameters
#mode is either 1 for mixture of point mass and single beta distn, or 2 for mixture of multiple betas
#return.est: returns only estimate if TRUE, and estimated p's and q when FALSE
binshrink=function(x,q,mode,return.est=FALSE){
  n=length(x)
  J=log2(n)
  titable=cxxParentTItable(x)
  tit=titable$TItable
  ptit=titable$parent
  pp=binash(tit,ptit,q,mode)$p
  est=reverse.pp(tit,pp,q,mode)
  if(return.est==TRUE){
    return(est)
  }else{
    return(list(est=est,p=pp,q=q))
  }
}
