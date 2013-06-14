### How to use this code
# Before anything else, install and load the package "mixfdr" from CRAN
# and source in "normal mixture.R"
# Suppose you have observations z
# Note that z is assumed to have sampling variance 1
# We want to get fdr estimates and effect size estimates
#
# First fit a mixture model to z ("mixFdr"). This has effect size
# estimates and fdrs for the z's. It also has parameter estimates for the model.

# If you need estimates for other points (like new data, or to plot)
# effectSize gives posterior mean, variance estimates
# fdrMixModel gives fdr estimates
# all of these work on some vector of values (new or old Z's)
# and need a fitted model for parameters

# more details are in the mixfdr package help.

# Example:
# We want to fit a mixture model with 
# J = 3 groups
# and the theoretical null

# the code to do this is
require(mixfdr)
 set.seed(100)
 z = c(rnorm(1000), rnorm(100,3,1)) # some example
 m = mixFdr(z, J = 3, theonull = TRUE)
 delta = m$effectSize
 varDelta = m$effectPostVar
 fdrs = m$fdr
# there are many more options - see the function mixmodel in "normal mixture" for details 
 
# Suppose we wanted fdr and effect size curves
# then we could do

 s = seq(-5,5,by=0.01)
 effCurve = effectSize(s,m)[,1]
 fdrCurve = fdrMixModel(s,m,abs(m$mu - m$mu[1])<=0.2) # this tells it which indices to consider null

# compare to the Bayes estimator for this problem
 mTrue = list(pi = c(10,1)/11, mu = c(0,3), sig = c(1,1))
# note that the model parameters are in terms of the MARGINAL not of the prior
 trueEff = effectSize(s,mTrue)[,1]
 trueFdr = fdrMixModel(s,mTrue,c(TRUE,FALSE))
 par(mfrow = c(2,1))
 plot(s, trueEff, t = 'l', main = "Effect Size - Black is true", xlab = "z", ylab = "E(delta|z)", ylim = c(-3,3))
 lines(s, effCurve, col = 2)
 plot(s, trueFdr, t = 'l', main = "fdr - Black is true", xlab = "z", ylab = "fdr", ylim = c(0,1))
 lines(s, fdrCurve, col = 2)
 
 ### simulation script

# the mixture model fitting code is fairly quick
# but each call to mixFdr tries a bunch of starts
# and the simulation is kind of large
# so this can still take a while to run...
#
# the paper's simulation results are saved in "oneSidedUnifEff.Rdata"
# and in the file "twoSidedUnifEff.Rdata""

## effect size
# source in simulationstuff.R
library(EbayesThresh)
resOneSided = simulation1full()
resTwoSided = simulation1full(twosided = TRUE)
# or load them in

# mean squared error tables
table1ebayess(resOneSided)
table1ebayess(resTwoSided)

# relative error assessment
library(abind)
restogether = abind(resOneSided, resTwoSided, along = 2)
table2ebayess(restogether)

# you can also get the individual inefficiency tables
# with for example
# table2ebayess(restogether[,1:4,,]) for the one sided case

#### Plotting code
library(RColorBrewer)
penalPalette3 = colorRampPalette(c("yellow","red"), space = "rgb")
a = brewer.pal(5, "Oranges")
penalPalette10 = colorRampPalette(c("cyan","blue"), space = "rgb") 
J3models = 9:13
J10models = 14:18

b1 = table1ebayess(restogether[,,,])
mat1 = cbind(b1[,,1], b1[,,2], b1[,,3])
mat1 = mat1[c(J3models,J10models,21),c(1:4,9:12,17:20,5:8,13:16,21:24)]
baseline = mat1[nrow(mat1),] # Bayes
mat1 = mat1 / matrix(baseline, nrow(mat1), ncol(mat1), byrow=T)
NAcol = rep(NA, nrow(mat1))
mat1 = cbind(mat1,NAcol)
mat1 = mat1[,c(1:4,25,5:8,25,9:12,25,13:16,25,17:20,25,21:24)]
x = c(1:4,NA,5:8,NA,9:12,NA,13:16,NA,17:20,NA,21:24)
matplot(x,t(mat1), t = 'l', lty = 1, col = c(penalPalette3(5), penalPalette10(5), brewer.pal(8,"Dark2")[8]), lwd = 1.5, main = "Parameter Choice Comparison", xlab = "", ylab = "Relative Error", xaxt="n", log = "y")
abline(v = 4.5, lty = 2)
abline(v = 8.5, lty = 2)
abline(v = 12.5, lwd = 2)
abline(v = 16.5, lty = 2)
abline(v = 20.5, lty = 2)
axis(1, at = 1:24, labels = rep(2:5,6), cex.axis = 0.75, padj = -1.5)
axis(3, at = 0.5+seq(2,22,by=4), tick = FALSE, labels = rep(c("Sparse (5)","Middle (50)","Dense (500)"),2), cex.axis = 0.75, padj = 2)
axis(1, at = c(6.5,18.5), labels = c("One sided","Two Sided"), padj = 1, tick = FALSE)
legend("topright", legend = c("J=3,no pen", "J=3,P=10","J=3,P=50","J=3,P=100","J=3,P=200","J=10,no pen", "J=10,P=10","J=10,P=50","J=10,P=100","J=10,P=200", "Bayes"), lty = 1, col = c(penalPalette3(5), penalPalette10(5), brewer.pal(8,"Dark2")[8]), lwd = 1.5,bg = "white")

# representatives of other methods compared to bayes
plotPal = brewer.pal(8, "Dark2")
models = c(16, 1, 3, 5, 7, 8, 19, 21)
names = c("Bayes", "MixFdr", "EBThresh", "SURE", "FDR","UnivSoft", "UnivHard","Spline")
b1 = table1ebayess(restogether[,,,])
mat1 = cbind(b1[,,1], b1[,,2], b1[,,3])
mat1 = mat1[models,c(1:4,9:12,17:20,5:8,13:16,21:24)]
baseline = mat1[8,] # Bayes
mat1 = mat1 / matrix(baseline, nrow(mat1), ncol(mat1), byrow=T)
NAcol = rep(NA, nrow(mat1))
mat1 = cbind(mat1,NAcol)
mat1 = mat1[,c(1:4,25,5:8,25,9:12,25,13:16,25,17:20,25,21:24)]
x = c(1:4,NA,5:8,NA,9:12,NA,13:16,NA,17:20,NA,21:24)



# subsetting 1:15 or 16:29 works to produce onesided or twosided plots 
# (subset the cols of mat1)
op = par(mar=c(5, 4, 4, 7)+.1, xpd = FALSE)
matplot(x[16:29],t(mat1[,16:29]), t = 'l', lty = 1, col = plotPal, main = "Comparison of Effect Size Estimation Methods", xlab = "", ylab = "Relative Error", xaxt="n", log = "y", ylim = c(1,6),lwd = 0.25+c(3,1,1,1,1,1,1,1))
abline(v = 4.5, lty = 2)
abline(v = 8.5, lty = 2)
abline(v = 12.5, lwd = 2)
abline(v = 16.5, lty = 2)
abline(v = 20.5, lty = 2)
axis(1, at = 1:24, labels = rep(2:5,6), cex.axis = 0.75, padj = -1.5)
axis(3, at = 0.5+seq(2,22,by=4), tick = FALSE, labels = rep(c("Sparse (5)","Middle (50)","Dense (500)"),2), cex.axis = 0.75, padj = 2)
axis(1, at = c(6.5,18.5), labels = c("One sided","Two Sided"), padj = 1, tick = FALSE)
par(xpd=TRUE)

# change to 24.44 for two sided
legend(24.44,6.95, legend = names, lty = 1, col = plotPal[c(8,1:7)],bg = "white",lwd = 0.25+c(1,2,1,1,1,1,1,1), bty = 'n')
par(op)