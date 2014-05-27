# FSR Example

This is an example of how the false sign rate and the false discovery rate can
be very different when a large proportion of the nulls are false.

In this example, we assume that the true $\beta_j$ ($j=1,\dots,n$) are normally distributed, $\beta_j \sim N(0,1)$. 
So no $\beta_j$ is exactly 0, and all the null hypotheses $H_j: \beta_j=0$ are false.
Since all the nulls are false, no discovery is a false discovery,
and the ``true" FDR is always 0, whatever threshold is used.



```r
opts_chunk$set(fig.path = "figure/lfsr/")
set.seed(100)
library(ashr)
```

```
## Loading required package: truncnorm
```

```r
library(qvalue)
n = 10000
beta = rnorm(n)
sebetahat = 1
betahat = beta + rnorm(n, 0, sebetahat)

z = betahat/sebetahat
p = pchisq(z^2, df = 1, lower.tail = FALSE)
z.ash = ash(z, 1, method = "fdr")
```

```
## [1] "normal likelihood"
```

```r
z.qvalue = qvalue(p)
```


In this example all the nulls are false, so the actual proportion of true nulls, $\pi_0$, is equal to 0. However, the estimated value of $\pi_0$ can be very different from zero. From qvalue we get
an estimated $\pi_0$ of 0.7146. This makes sense when you look at the
distribution of p values, and recall that qvalue makes the assumption that all $p$ values near 1 are null (see histogram below). In contrast, ash produces a much smaller estimate of $\pi_0$ 
($\pi0=$ 0.1258), because it allows that many $p$ values near 1 may actually arise
even when $\beta_j \neq 0$. Of course, in this case both qvalue and ash provide "conservative"
estimates of $\pi_0$, but the qvalue estimate is much more conservative. 



```r
hist(p, prob = TRUE)
abline(h = z.qvalue$pi0, lwd = 2)
text(0.95, z.qvalue$pi0 + 0.06, "pi0 (q value)")
abline(h = get_pi0(z.ash), col = 2, lwd = 2)
text(0.95, get_pi0(z.ash) + 0.06, "pi0 (ash)", col = 2)
```

![plot of chunk unnamed-chunk-2](figure/lfsr/unnamed-chunk-2.png) 


Because of its very conservative estimate of $\pi_0$, qvalue also
provides much more conservative estimates of the FDR. For example, at a $p$ value threshold of 0.05, qvalue estimates the FDR to be 0.2138 whereas ash estimates the FDR
as 0.0405. Of course, the true value of the FDR is zero, so both are conservative. 

Even though the true FDR is zero at any threshold, it seems, intuitively, that "discoveries" with
a $p$ value near 1 are not very meaningful. One way to formalize this is to consider the ``false sign rate" instead of the false discovery rate. The local false sign rate (lfsr) for each discovery is the probability that the true $\beta_j$ has the opposite sign to the corresponding $z$ score, $z_j$. The tail false sign rate (FSR) at a given threshold is the average of the lfsr values for discoveries that exceed that threshold. In this example, discoveries with a $p$ value near 1 will have a high lfsr, even though they have a local false discovery rate (lfdr) of 0. When, as here, $\pi_0$ is large, the lfsr is much more meaningful than
the lfdr. When $\pi_0$ is small both measures are meaningful, 
and indeed are generally quite similar. Consequently we recommend lfsr for general usage.

In this example, at a $p$ value threshold of 0.05, $\ash$ estimates the local false sign rate (lfsr) to be 0.1471 and the tail false sign rate (FSR) to be 0.0786. The actual FSR at that threshold is 0.0557.



