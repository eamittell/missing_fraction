### Simulations of true values for selection across episodes for the missing fraction problem
# Lizy

library(mvtnorm)

betas <- c(-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.05,0,0.05,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
cvrs <- c(0,0.2,0.4)
S_logistic <- function(a,b,m,s){
	low <- m - 6*s; up <- m + 6*s; # get the sample range
	fun <- function(z,a,b,m,s){(exp(a+b*z)/(1+exp(a+b*z)))*dnorm(z,m,s)} # function to get mean absolute fitness
	muW <- integrate(f=fun,lower=low,upper=up,a=a,b=b,m=m,s=s)$value # obtain mean absolute fitness
	ffun <- function(z,a,b,m,s){(z-m)*(exp(a+b*z)/(1+exp(a+b*z)))*dnorm(z,m,s)} # function to get covariance between trait and absolute fitness
	return((integrate(f=ffun,lower=low,upper=up,a=a,b=b,m=m,s=s)$value-m*muW)/muW) # normalise to relative fitness to get the covariance between trait and relative fitness, i.e. the selection differential
}

# Object to save simulated true values from each combination of parameters
true <- NULL
for(i in 1:length(betas)){
# Which beta using in this simulation
beta1 <- betas[i]
# Object to store each specific simulated result
temp2 <- NULL
# Run for range of phenotypic covariances for each viability selection regime
for(k in 1:length(cvrs)){
temp <- NULL
# Which phenotypic covariance using in this simulation
cvr <- cvrs[k]

# P matrix in the absence of selection
V1 <- 1
V2 <- 1
rho <- 0.4 + cvr
P <- matrix(c(V1,rho*sqrt(V1)*sqrt(V2),rho*sqrt(V1)*sqrt(V2),V2),2,2)

# Viability selection on early-life trait
S1a_log <- S_logistic(0,beta1,0,1)

# Prior viability selection on later-life trait
S2a_log <- S1a_log/(P[2,2])*P[1,2]

# Variance of survivors from first to second episode
V2b <- V2 - S2a_log^2

# Fecundity selection on later-life trait of survivors
beta2 <- 0.25
S2b <- beta2*V2b

temp <- c(beta1,beta2,cvr,V1,V2,rho,S1a_log,S2a_log,V2b,S2b)
temp2 <- rbind(temp2,temp)
}
true <- rbind(true,temp2)
}
