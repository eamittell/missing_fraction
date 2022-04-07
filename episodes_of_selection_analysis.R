# Composite function to calculate the prior viability selection differential from the 
# logistic intercept and slope (a and b), trait mean and sd (m and s)
S_logistic <- function(a,b,m,s){
# get the sample range
	low <- m - 6*s; up <- m + 6*s;
# function to get absolute fitness
	fun <- function(z,a,b,m,s){(exp(a+b*z)/(1+exp(a+b*z)))*dnorm(z,m,s)}
# obtain absolute fitness
	muW <- integrate(f=fun,lower=low,upper=up,a=a,b=b,m=m,s=s)$value
# function to get covariance between absolute fitness and trait 2
	ffun <- function(z,a,b,m,s){(z-m)*(exp(a+b*z)/(1+exp(a+b*z)))*dnorm(z,m,s)}
# normalise to relative fitness
	return((integrate(f=ffun,lower=low,upper=up,a=a,b=b,m=m,s=s)$value-m*muW)/muW)
}

# Infinite population
# Viability selection on trait 1
S1a_inf <- S_logistic(0,beta,0,sqrt(P[1,1]))
# calculation of what prior viability selection on trait 2 should be
S2a_inf <- S1a_inf/(P[1,1])*P[1,2]

# Simulated population
# Phenotypic S on trait 1
S1a <- weighted.mean(d$z1,d$W1)-mean(d$z1)
# Phenotypic prior sel diff trait 2
S2a <- weighted.mean(d$z2,d$W1)-mean(d$z2)
# total selection of trait 2, but with missing fraction problem -- i.e. fecundity selection on z2
d6 <- subset(d,d$W1==1)
S2b <- weighted.mean(d6$z2,d6$W)-mean(d6$z2)
# total selection of trait 2, if we had perfect knowledge
S2t <- weighted.mean(d$z2,d$W)-mean(d$z2)

library(norm2)
# The variances and covariances estimated using expectation maximum likelihood
ML_P <- emNorm(cbind(d$z1,d$z2_obs))
sigma2 <- summary(ML_P)$param$sigma
# Covariance-variance matrix z1 and values observed for z2 
sigma2

# True
covar <- cov(d[,c("z1","z2","z2_obs")],use="pairwise.complete.obs")

### Same analyses and results for half-sibs ###
