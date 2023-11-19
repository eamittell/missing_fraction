# Covariance-variance matrix z1, all z2 and values observed for z2
covar <- cov(d[,c("z1","z2","z2_obs")],use="pairwise.complete.obs")

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

# Phenotypic viability selection on trait 1
# In an infinite population across known parameters
S1a_inf <- S_logistic(0,beta,0,sqrt(P[1,1]))
# Calculation of what prior viability selection on trait 2 should be
S2a_inf <- S1a_inf/(P[2,2])*P[1,2]

# Phenotypic S on trait 1 in simulated data
# Known if measured early-life traits and viability selection in a focal population
S1a <- weighted.mean(d$z1,d$W1)-mean(d$z1)

library(norm2)
# The variances and covariances estimated using expectation maximum likelihood
ML_P <- emNorm(cbind(d$z1,d$z2_obs))
sigma2 <- summary(ML_P)$param$sigma
# Covariance-variance matrix z1 and values observed for z2 
sigma2

# Estimated phenotypic prior selection differential trait 2
S2a_est <- (S1a/sigma2[2,2])*sigma2[1,2]

# Variance of survivors from first to second episode
V2b <- P[2,2] - S2a_logˆ2

# Fecundity selection differential on later-life trait of survivors
# Remembering that selection differential = selection gradient (beta)*variance
S2b <- beta2*V2b
# Total selection of trait 2
S2t <- S2a_log + S2b
# Total estimated selection of trait 2
S2t_est <- S2a_est + S2b
