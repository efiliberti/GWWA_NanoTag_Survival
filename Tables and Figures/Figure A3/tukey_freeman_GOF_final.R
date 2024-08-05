
## Figure A2 in manuscript. This script runs a Tukey-Freeman goodness-of-fit 
## test on our simplest model.  

## Required packages
library(IPMbook)
library(dplyr)
library(tidyr)
library(ggplot2)
library(EnvStats)
library(sysfonts)
library(nimble)
library(MCMCvis)
library(gridExtra)
library(jagsUI)

################################################################################
####                         GOODNESS-OF-FIT TEST                           ####
################################################################################

## Build the model
sink("cjs.tag.sex.eff.jags")
cat("model{
  # Multinomial likelihood
  for (t in 1:(n.occasions-1)){  # rows of the m-array/p-array
    marr[t, ] ~ dmulti(parr[t, ], R[t])
  }

  # Priors
  # Using vectors for phi and p simplies creating the p-array
  for (t in 1:(n.occasions - 1)){
    phi[t] <- phi0
    p[t] <- p0
  }
  phi0 ~ dbeta(1, 1)         # Prior for survival
  p0 ~ dbeta(1, 1)           # Prior for recapture

  # Create the p-array, the cell probabilities of the m-array:
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    parr[t,t] <- phi[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      parr[t,j] <- prod(phi[t:j]) * prod(1 - p[t:(j-1)]) * p[j]
    }
    # Below main diagonal
    for (j in 1:(t-1)){
      parr[t,j] <- 0
    }
    # Last column: probability never recaptured
    parr[t,n.occasions] <- 1 - sum(parr[t,1:(n.occasions-1)])
  }


# GoF for cap.-recap. data: Freeman-Tukey test statistics 
for (t in 1:(n.occasions-1)){
# Simulated m-arrays
marr.pred[t,1:n.occasions] ~ dmulti(parr[t,], R[t])

# Expected values and test statistics
for (j in 1:n.occasions){
marr.E[t,j] <- parr[t,j] * R[t]
E.org[t,j] <- pow((pow(marr[t,j], 0.5) - pow(marr.E[t,j], 0.5)), 2)
E.new[t,j] <- pow((pow(marr.pred[t,j], 0.5) - pow(marr.E[t,j], 0.5)), 2)
  } #j
} #t

fit <- sum(E.org) 
fit.new <- sum(E.new) 

}
", fill = T)


sink()


## Bundle data
known.state.cjs	<- function(ch){
  state	<- ch
  for	(i	in	1:dim(ch)[1]){
    n1	<- min(which(ch[i,]==1))
    n2	<- max(which(ch[i,]==1))
    state[i,n1:n2]	<- 1
    state[i,n1]	<- NA
  }
  state[state==0]	<- NA
  return(state)
}

marr<-marray(y)
R <- rowSums(marr)
jags.data <- list(marr = marr, n.occasions=ncol(marr), R = R)


known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    #state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}


## Initial Values
inits <- function(){list()}

## Parameters Monitored
parameters <- c("phi",  "p",   "fit", "fit.new")

## MCMC settings
ni <- 100000
nt <- 2
nb <- 20000
nc <- 3

## Run and summarize the model

## call jags from R
cjs <- jags(jags.data, inits, parameters, "cjs.tag.sex.eff.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


## Summarize Posteriors
print(cjs, digits = 3)

# Evaluation of Fit
plot(cjs$sims.list$fit, cjs$sims.list$fit.new, xlab = "Discrepancy
actual data", ylab = "Discrepancy replicate data", las = 1,
ylim = c(0, 10), xlim = c(0, 10), bty ="n")
abline(0, 1, col = "black", lwd = 2)
mean(cjs$sims.list$fit.new > cjs$sims.list$fit)
