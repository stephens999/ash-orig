The following R code simulates some data where the effects are normally distributed, and the first 500 observations have good precision, and the next 500 have poor precision.


```r
# install q value package source('http://bioconductor.org/biocLite.R')
# biocLite('qvalue')
library("qvalue")
library("lattice")  #library for some of the plots

# set up some data with mixture of values of s
set.seed(100)
s.good = 0.5
s.poor = 10
J.good = 500
J.poor = 500
J = J.good + J.poor
beta = c(rnorm(J, 0, 1))
s = c(rep(s.good, J.good), rep(s.poor, J.poor))
betahat = beta + rnorm(J, 0, s)
# compute the usual zscore and corresponding p value
zscore = betahat/s
pval = pchisq(zscore^2, df = 1, lower.tail = F)
```

As expected the $p$ values from the poor observations are approximately uniform, whereas those from the good observations  are enriched for small $p$ values:


```r
p = data.frame(pval = pval, type = c(rep("GOOD", 500), rep("POOR", 500)))
require(ggplot2)
```

```
## Loading required package: ggplot2
## 
## Attaching package: 'ggplot2'
## 
## The following object is masked from 'package:qvalue':
## 
##     qplot
```

```r
pdf("figures/good_poor_hist.pdf", height = 3, width = 6)
ggplot(p, aes(x = pval)) + facet_grid(. ~ type) + geom_histogram(aes(y = ..density..), 
    binwidth = 0.05, colour = "black", fill = "white")

# old version, using trellis
# histogram(~pval|type,p,layout=c(1,2),xlim=c(0,1),breaks=seq(from=0,to=1,length=20),type='count')
dev.off()
```

```
## pdf 
##   2
```


Now what happens if we apply FDR methods to the all the data is
that the uniform $p$ values from the poor observations
add noise relative to looking at the good observations only. This impacts both the estimate of $pi_0$ from qvalue, and the number of findings that are significant at a given FDR.

```r
qq.all = qvalue(p$pval)
qq.good = qvalue(p$pval[p$type == "GOOD"])
print(c(qq.good$pi0, qq.all$pi0))
```

```
## [1] 0.3368 0.6932
```

```r
print(c(sum(qq.good$qvalues < 0.05), sum(qq.all$qvalues < 0.05)))
```

```
## [1] 193 114
```


In contrast, if you use adaptive shrinkage, it makes no difference whether you use all the observations or the just the good ones -
the noise is ignored (as it should be!)


```r
library("ashr")
```

```
## Loading required package: truncnorm
## Loading required package: SQUAREM
```

```r
library("locfdr")
```

```
## Loading required package: splines
```

```r
beta.ash.all = ash(betahat,s,method="fdr")
```

```
## [1] "normal likelihood"
```

```r
beta.ash.good = ash(betahat[p$type=="GOOD"],s[p$type=="GOOD"],method="fdr")
```

```
## [1] "normal likelihood"
```

```r
beta.ash.all.z = ash(betahat/s,1,method="fdr")
```

```
## [1] "normal likelihood"
```

```r
beta.ash.good.z = ash(betahat[p$type=="GOOD"]/s[p$type=="GOOD"],1,method="fdr")
```

```
## [1] "normal likelihood"
```

```r
beta.locfdr.all.z = locfdr(betahat/s,nulltype=0,plot=0)
beta.locfdr.good.z = locfdr(betahat[p$type=="GOOD"]/s[p$type=="GOOD"],nulltype=0,plot=0)
res = rbind(data.frame(x=beta.ash.good$lfsr,y=beta.ash.all$lfsr[p$type=="GOOD"],type="ash"), 
            #data.frame(x=beta.ash.good.z$lfsr,y=beta.ash.all.z$lfsr[p$type=="GOOD"],type='ash, z scores'),
            data.frame(x=beta.locfdr.good.z$fdr,y=beta.locfdr.all.z$fdr[p$type=="GOOD"],type='locfdr'), 
            data.frame(x= qq.good$qvalues, y = qq.all$qvalues[p$type=="GOOD"],type="qvalue") )
  
pdf("figures/good_vs_all.pdf",height=3,width=6.5)
  pp= ggplot(data=res,aes(x,y)) +geom_point(shape=1) +
        facet_grid(. ~ type) +
      geom_abline(colour = "black") +
        xlab("lfsr/lfdr/qvalue (Good data only)") +
        ylab("lfsr/lfdr/qvalue (All data)")


  print(pp +scale_y_continuous(limits=c(0,1)) +
          scale_x_continuous(limits=c(0,1))  +
          coord_equal(ratio=1))
  
dev.off()
```

```
## pdf 
##   2
```


