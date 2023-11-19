# Libraries required
library(mvtnorm)

## SIMULATION SETUP HALF SIBS ##
## Half-sib data are used in the manuscript
## However, we have carried out the same analyses on full-sib data and found no difference using the two different relatedness structures
## The code for full-sib data generation are included in this repository

##### Variables to use if running full loops
# Residual covariance -> phenotypic covariance trait 2 for example 0.4(G) + 0.4(R) = 0.8(P)
#cvrs <- c(0,0.4) -- Put range of residual (and so phenotypic) covariances that you want to simulate data for
# Range of strengths of selection on early-life trait 1
#betas <- c(-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.05,0,0.05,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# For one example simulation
cvr <- 0.4
beta <- 0.9
n.sires <- 500 
n.offspring <- 5

# Generate phenotypes
d <- expand.grid(sire=1:n.sires,o=1:n.offspring)
d <- d[order(d$sire),]
d$sire <- paste("S",d$sire,sep="")
d$dam <- paste("D",1:(n.sires*n.offspring),sep="") # Half sibs
d$id <- paste("O",1:(n.sires*n.offspring),sep="")

# Half sibs 
n.dams <- n.sires*n.offspring

# Organise the data into a pedigree, add dam and sire records before offspring
pedigree <- d[,c("id","sire","dam")]
dam.p <- data.frame(id=unique(d$dam),sire=NA,dam=NA)
pedigree <- rbind(dam.p,pedigree)
sire.p <- data.frame(id=unique(d$sire),sire=NA,dam=NA)
pedigree <- rbind(sire.p,pedigree)

# Genetic covariance 0.4
G <- matrix(c(0.5,0.4,0.4,0.5),2,2)

# Residual variance-covariance matrix
R <- matrix(c(0.5,cvr,cvr,0.5),2,2)

# Phenotypic var-covar matrix
P <- G + R

# Breeding values for each trait
a <- rbv(pedigree, G)
a1 <- a[3001:5500,] # Half sibs
d$a1 <- a1[,1]
d$a2 <- a1[,2]

# Check 1
a[3001:3010,] # Half sibs
d[1:10,]

# Check 2: offspring on single parent regression of a should be 0.5 (so roughly 1 if we multiple by 2)
d$father_a1 <- a[match(d$sire,dimnames(a)[[1]]),1]
d$mother_a1 <- a[match(d$dam,dimnames(a)[[1]]),1]
summary(lm(a1~father_a1,data=d))$coefficients[2,]
coef(lm(a1~father_a1,data=d))[2]*2
summary(lm(a1~mother_a1,data=d))$coefficients[2,]
coef(lm(a1~mother_a1,data=d))[2]*2

# Residual effects
r <- rmvnorm(n.sires*n.offspring,c(0,0),R)
rownames(r) <- d$id # Add id
d$r1 <- r[,1]
d$r2 <- r[,2]

# Compose phenotypes
d$z1 <- d$a1 + d$r1
d$z2 <- d$a2 + d$r2

# Survival based on early-life trait 1
inv.logit <- function(x){exp(x)/(1+exp(x))}
surv_func <- function(x,a=0,b=beta){inv.logit(a+b*x)}
d$W1 <- rbinom(dim(d)[1],1,surv_func(d$z1)) # Fitness at viability stage

# Later-life trait phenotypes
d$z2_obs <- d$z2
d$z2_obs[which(d$W1==0)] <- NA # Don't observe later-life trait for those individuals that do not survive viability stage

# Later-life trait breeding values
d$a2_obs <- d$a2
d$a2_obs[which(d$W1==0)] <- NA

# Fitness consequence for trait 2
d$W2 <- rpois(length(d$z1),exp(0+0.25*d$z2))
# Total lifetime fitness
d$W <- d$W2*d$W1
# Fitness observed for trait 2
d$W_obs <- d$W

str(d)
summary(d)
