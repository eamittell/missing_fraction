# Load libraries
library(rjags)

# Results were the same from a probit function instead of the logistic which is shown

# Jags simulation model -- genetic logistic
missingModgen <- "model{
# Priors
# Genetic
sd_a ~ dunif(0,3)
tau_a <- 1/sd_a^2
# Residual
res_sd ~ dunif(0,3)
res_tau <- 1/res_sd^2

mu ~ dnorm(0,0.0001)
alpha ~ dnorm(0,0.0001)
beta ~ dnorm(0,0.0001)

# Parents breeding values
for(j in 1:n_s){
	a_s[j] ~ dnorm(0,tau_a)
}
for(j in 1:n_d){
	a_d[j] ~ dnorm(0,tau_a)
}

# Simulation
for(i in 1:n_ind){
	# Genotypes
	a[i] ~ dnorm((a_d[dam[i]]+a_s[sire[i]])/2,tau_a*2)
	# Phenotypes
	z[i] ~ dnorm(mu + a[i],res_tau)
	# Observed phenotypes with little error
	z_obs[i] ~ dnorm(z[i], 1000)
	# Survival probability based on the fitness function
	s[i] ~ dbern(ilogit(alpha + beta*a[i]))
}

# Track variances
V_a <- sd_a^2
V_res <- res_sd^2
Vp0 <- sd_a^2 + res_sd^2
}"
writeLines(missingModgen,"./missingModgen.jags")

dat2 <- list(
n_ind=dim(d)[1],
n_d=n.dams,
n_s=n.sires,
z_obs=d$z2_obs,
dam=as.integer(as.factor(d$dam)),
sire=as.integer(as.factor(d$sire)),
s=d$W1
)

m2 <- jags.model(file="./missingModgen.jags",data=dat2,n.chains=2,n.adapt=1000,quiet=TRUE)
s2 <- jags.samples(model=m2,variable.names=c("mu","alpha","beta","V_a","V_res","Vp0","a"),n.iter=200000,thin=100,quiet=TRUE)

# Function to estimate prior viability selection
S_logistic <- function(a,b,m,s){
	low <- m - 6*s; up <- m + 6*s; # get the sample range
	fun <- function(z,a,b,m,s){(exp(a+b*z)/(1+exp(a+b*z)))*dnorm(z,m,s)} # function to get fitness
	muW <- integrate(f=fun,lower=low,upper=up,a=a,b=b,m=m,s=s)$value # obtain fitness
	ffun <- function(z,a,b,m,s){(z-m)*(exp(a+b*z)/(1+exp(a+b*z)))*dnorm(z,m,s)} # function to get relative fitness
	return((integrate(f=ffun,lower=low,upper=up,a=a,b=b,m=m,s=s)$value-m*muW)/muW) # normalise to relative fitness
}

betas <- c(-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.05,0,0.05,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# Infinite population genetic viability selection
# Viability selection on trait 1
S1ag_log <- S_logistic(0,beta,0,1)
# calculation of what prior viability selection on trait 2 should be in infinite population
S2ag_log <- (G%*%c(S1ag_log/P[1,1],0))[2]

# In each finite population run
# Genetic viability selection on trait 2
post.S2ag_jags <- array(dim=2000)
for(i in 1:2000){
	post.S2ag_jags[i] <- S_logistic(s2$alpha[1,i,1],s2$beta[1,i,1],mean(s2$a[1:2500,i,1]),sqrt(s2$V_a[1,i,1]))
}
S2ag_jags <- mean(post.S2ag_jags)

# Later-life trait (2)
# True viability induced change in breeding value
S2a_gen <- weighted.mean(d$a2,d$W1)-mean(d$a2)
# True fecundity induced change in breeding value
S2b_gen <- weighted.mean(d$a2,d$W2)-mean(d$a2)
# True total change in breeding value
S2t_gen <- weighted.mean(d$a2,d$W)-mean(d$a2)
