# This code performs the analysis in Section 2 of the paper
# mu = arcsin(sqrt(delta)); delta = prob

# Here is an explanation of the beta-binomial model code
# "emfitter" fits a J component beta prior to data z
# throughout, z are observations and NN is a vector with the number of trials for each z
# that is, z[i] ~ Binom(delta[i], NN[i])
# so we have variable names like ztr(aining), NNtr(aining) etc

# the binomial stuff is rudimentary
# emfitter fits the model. it doesn't have many parameters beyond the number of groups J.
# Unlike the normal mixture case, penalization and constraints are not implemented.

# emfitter returns a model object that other function use
# binomialEffectSize gives posterior means (for delta)
# binomMixDens gives the marginal density
# binomgroupProbs gives posterior group probabilities
# EmuZ gives E(mu|z)

# baseball data
data = read.csv("basedat.txt")
# binomial functions
source("binomialFuncs.R")

pitcher = as.vector(data[,3])
train.ab = as.matrix(data[,5:7])
train.h = as.matrix(data[,11:13])
test.ab = as.matrix(data[,8:10])
test.h = as.matrix(data[,14:16])

NNtr = rowSums(train.ab)
ztr = rowSums(train.h)
NNte = rowSums(test.ab)
zte = rowSums(test.h)

indtr = (NNtr >= 11) 
indte = indtr & (NNte >= 11)

estJS = function(ztr,NNtr){
	X = asin(sqrt((ztr+.25)/(NNtr+.5)))
	sigs = .25/NNtr
	muHat = sum(X/sigs)/sum(1/sigs)
	shrink = 1 - (length(ztr) - 3)/(sum((X-muHat)^2/sigs))
	shrink = max(0,shrink)
	return(muHat + shrink*(X-muHat))	
}

predictionError = function(fits, zte,NNte,ztr,NNtr){
	target = asin(sqrt((zte+0.25)/(NNte+0.5)))
	SSPE = sum((fits-target)^2) - 0.25 * sum(1/NNte)
	nfits = asin(sqrt((ztr+0.25)/(NNtr+0.5)))
	SSNa = sum((nfits-target)^2) - 0.25 * sum(1/NNte)
	return(SSPE/SSNa)
}

e = emfitter(ztr[indtr], NNtr[indtr], J = 3)
fits = EmuZ(ztr[indte], NNtr[indte], e)
predictionError(fits, zte[indte], NNte[indte], ztr[indte], NNtr[indte])
#error for single overall model

#split up into pitchers and nonpitchers
ind1 = pitcher==1 & indtr
ind2 = pitcher==0 & indtr
e1 = emfitter(ztr[ind1], NNtr[ind1], J = 3, maxiter = 1000)
e2 = emfitter(ztr[ind2], NNtr[ind2], J = 3, maxiter = 1000)
newFits = ztr
newFits[ind1  & indte] = EmuZ(ztr[ind1 & indte],NNtr[ind1 & indte],e1)
newFits[ind2  & indte] = EmuZ(ztr[ind2 & indte],NNtr[ind2 & indte],e2)
predictionError(newFits[indte], zte[indte], NNte[indte], ztr[indte], NNtr[indte]) #overall fit of two models
predictionError(newFits[ind1 & indte], zte[ind1 & indte], NNte[ind1  & indte], ztr[ind1 & indte], NNtr[ind1 & indte]) #pitchers
predictionError(newFits[ind2 & indte], zte[ind2 & indte], NNte[ind2 & indte], ztr[ind2 & indte], NNtr[ind2 & indte]) #others

# try to find the sd of the performance difference under simulation
# generate new test data one and for all
# then try new training data sets
# all with the estimated batting averages
# there's some mucking around with index sets to make it easy to do pitchers, batters, etc
# compare to the james-stein estimator
# or the overall mean
set.seed(500)
iset = ind1 & indte
eiset = e1
testData = rbinom(sum(iset), NNte[iset], binomialEffectSize(ztr[iset],eiset,NNtr[iset]))
predErrDiff = 1:100
for(i in 1:100){
	ztrNew = rbinom(sum(iset), NNtr[iset], binomialEffectSize(ztr[iset],eiset,NNtr[iset]))
	eNew = emfitter(ztrNew,NNtr[iset], J = 3, maxiter = 1000)
#	Xnew = asin(sqrt((ztrNew+.25)/(NNtr[iset]+.5)))
	p1 = predictionError(EmuZ(ztrNew,NNtr[iset],eNew),testData,NNte[iset],ztrNew,NNtr[iset])
	p2 = predictionError(estJS(ztrNew,NNtr[iset]),testData,NNte[iset],ztrNew,NNtr[iset])
	predErrDiff[i] = p1 - p2	
}
sd(predErrDiff)
# results are similar for overall mean
# normalized performance (and sd) vary a lot by test data, not so much by training data

# difference on batters is probably stable, difference on pitchers is probably not

#now add ash to comparison
require(ashr)
phat.tr=(1+ztr)/(2+NNtr)
se.tr = sqrt(phat.tr*(1-phat.tr)/(NNtr+1))
phat.ash = ash(phat.tr-median(phat.tr),se.tr)
plot(phat.tr,phat.ash$PosteriorMean+median(phat.tr))
abline(a=0,b=1)

estAsh = function(ztr,NNtr){
  X = asin(sqrt((ztr+.25)/(NNtr+.5)))
  sigs = .25/NNtr
  muHat = sum(X/sigs)/sum(1/sigs)
  return(ash(X-muHat,sqrt(sigs))$PosteriorMean+muHat)
}

indtr=indte
xx.ash=estAsh(ztr[indtr],NNtr[indtr])

xx.JS = estJS(ztr[indtr],NNtr[indtr])
plot(x[indtr],xx.ash)


predictionError(xx.ash, zte[indtr],NNte[indtr],ztr[indtr],NNtr[indtr])
predictionError(xx.JS, zte[indtr],NNte[indtr],ztr[indtr],NNtr[indtr])

x.te= asin(sqrt((zte+0.25)/(NNte+0.5)))
mean((x.te[indtr]-xx.ash)^2)
mean((x.te[indtr]-xx.JS)^2)