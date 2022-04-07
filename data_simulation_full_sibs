# Libraries required
library(mvtnorm)

## SIMULATION SETUP FULL SIBS ##
n.sibs <- 5 # Number of sibs per sibship (offspring of each parent)
n.sires <- 500 # Number of families (i.e. sample size, which could be number of occupied nest-boxes)
alpha <- 0 # Logistic regression intercept for z1 (early-life trait)
beta <- 0.9 # Controls the strength of viability selection on z1
# A range of beta were used in simulations: 
betas <- c(-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.05,0,0.05,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
# Genetic covariance set at 0.4
S0 <- matrix(c(0.5,0.4,0.4,0.5),2,2) # Genetic var-covar matrix between z1 and z2 (later-in-life trait)
# Phenotypic covariance set at 0.8 for highly correlated traits, 
# 0.4 for more weakly correlated traits
P0 <- matrix(c(1,0.8,0.8,1),2,2) # Phenotypic var-covar matrix

# Generate phenotypes
d <- expand.grid(sire=1:n.sires,sib=1:n.sibs)
s <- rmvnorm(n.sires,c(0,0),S0) # Genetic component of the two traits for each family
z <- s[d$sire,] + rmvnorm(n.sires*n.sibs,c(0,0),P0-S0) # Phenotypes for each individual 
# in the offspring generation
d$z1 <- z[,1] # Early-life trait
d$z2 <- z[,2] # Later-life trait

# Do individuals survive to express z2 or not?
inv.logit <- function(x){exp(x)/(1+exp(x))} # Define the inverse logit
FF1 <- function(z){rbinom(length(z),1,inv.logit(alpha+beta*z))} # Prior survival function for z1
# Generates binomial output for surviving (1) or dying (0) from early-life viability selection
# As alpha is set at 0 then 50% of the population will survive/die

# Prior survival from viability selection on z1
d$W1 <- FF1(z[,1])
# Consequence for z2
d$W2 <- rpois(length(z[,1]),exp(0+0.25*z[,2])) 
# Poisson generates a range of numbers of offspring, and therefore, fitness at the later-life stage

# Total lifetime fitness
d$W <- d$W2*d$W1
# Fitness observed for z2
d$W_obs <- d$W2
# z2 has no fitness when the individual died in the first episode
d$W_obs[which(d$W1==0)] <- NA

# The phenotypes as observed in the focal episode
# (i.e. phenotypes of those already dead must be NA)
d$z2_obs <- d$z2
d$z2_obs[which(d$W1==0)] <- NA
