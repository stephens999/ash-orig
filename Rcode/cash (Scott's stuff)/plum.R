### glm p.179
#data <- c(98, 0, 0, 51, 2, 1, 34, 6, 3, 35, 5, 8,
#	32, 10, 9, 23, 7, 8, 12, 6, 10, 4, 2,5)
#yrs <- c(5.8, 15.0, 21.5, 27.5, 33.5, 39.5, 46.0, 51.5)
#n <- 8
#k <- 3
#data <- matrix(data, n, k, byrow=T)
#fit1 <- plum(data, log(yrs))
#fit2 <- plum(data, cbind(yrs,log(yrs)))
#fit3 <- plum(data, log(yrs), link="c-log log")
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#nrows <- 4
#ncols <- 9
#freq <- scan("http://www.stat.uchicago.edu/~pmcc/glm/GLM175.txt")
#data <- matrix(freq, ncol=ncols, byrow=T)
#cheese <- as.factor(1:4)
#fit0 <- plum(data)
#fit1 <- plum(data, cheese)
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plum <- function(data, X, link="logit", maxcyc=15, tol=1.0e-4){
#
# fits the proportional odds model or related model
# the matrix X should be of full rank, as should the augmented matrix 1X
#
if(is.factor(data)) data <- model.matrix(~data-1)
n <- dim(data)[1]
k <- dim(data)[2]
nk <- n*k
nk1 <- nk - n
p <- 0
if(missing(X)) X <- matrix(rep(1, n), n, 1)
if(!missing(X)){
	if(is.factor(X)) X <- model.matrix(~X)[, 2:length(levels(X))]
	if(is.vector(X)) X <- matrix(X, n, 1)
	}
qr <- qr(cbind(rep(1, n), X))
X <- X[, qr$pivot[2:qr$rank] - 1]
p <- qr$rank - 1 
m <- data %*% rep(1, k)
tot <- t(rep(1, n)) %*% data
ctot <- cumsum(tot)[1:(k-1)]/sum(m)
# initial values
beta <- c(log(ctot/(1-ctot)), rep(0, p))
if(link == "probit") beta <- c(qnorm(ctot), rep(0, p))
if(link == "c-log log") beta <- c(log(-log(1-ctot)), rep(0, p))
if(link == "log log") beta <- c(-log(-log(ctot)), rep(0, p))
fitted <- m %*% tot / sum(tot)
dev0 <- 2*sum(data * log(pmax(data, 1)/fitted))

theta <- gl(k-1, 1, nk1)
X0   <- model.matrix(~ theta - 1)
dupl <- model.matrix(~ gl(n, k-1, nk1) - 1)
if(p > 0) X <- cbind(X0, dupl %*% X) else X <- X0

mk <- rep(m, rep(k, n))
m <- rep(m, rep(k-1, n))
y <- t(apply(data, 1, cumsum))[, 1:(k-1)]
y <- as.vector(t(y))
diff <- diag(k)
for(i in 2:k) diff[i-1, i] <- -1
ndg <- as.logical(rep(c(rep(1, k-1), 0), n))

for(cycle in 1:maxcyc){
	eta <- as.vector(X %*% beta)
	gamma <- 1/(1 + exp(-eta))
	if(link == "probit")    gamma <- pnorm(eta)
	if(link == "c-log log") gamma <- 1 - exp(-exp(eta))
	if(link == "log log")   gamma <- exp(-exp(-eta))
	cfit <- m * gamma
	mgamma <- cbind(matrix(gamma, n, k-1, byrow=T), rep(1, n))
	pi <- as.vector(t(mgamma %*% diff))
	fitted <- matrix(mk * pi, n, k, byrow=T)
	dev <- 2*sum(data * log(pmax(data, 1)/fitted))
	W <- diag(1.0/pi + 1/pi[c(2:nk, 1)])
	for(i in 2:nk){
		W[i-1, i] <- -1/pi[i]
		W[i, i-1] <- -1/pi[i]
		}
	der <- gamma*(1-gamma)
	if(link == "probit") der <- dnorm(eta)
	if(link == "c-log log") der <- exp(eta - exp(eta))
	if(link == "log log") der <- exp(-eta - exp(-eta))
	W <- diag(der) %*% W[ndg, ndg] %*% diag(der)
	deriv <- t(X) %*% W %*% ((y - cfit)/der)
	FI <- t(X) %*% diag(m) %*% W %*% X
	beta <- beta + solve(FI, deriv)
	if(cycle > 1 && abs(dev0 - dev) < tol) break else dev0 <- dev
	}	
list(coef=beta, fitted=fitted, deviance=dev, beta.cov=solve(FI), cycle=cycle,
	df=nk1-qr(X)$rank)
}