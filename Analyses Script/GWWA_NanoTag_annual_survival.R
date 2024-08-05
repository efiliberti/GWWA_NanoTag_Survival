## This script runs apparent annual survival CJS models on Golden-winged Warbler 
## CMR datasets collected between 2021-2023. Since site as a fixed effect on
## detection probability explained the most variation, survival models are run 
## with both fixed detection probability and site as a fixed effect on detection
## probability. A WAIC model comparison is then performed to discern which 
## covariates explain the most variation in survival rates. 


## Load packages
library(IPMbook)
library(dplyr)
library(tidyr)
library(ggplot2)
library(EnvStats)
library(sysfonts)
library(nimble)
library(MCMCvis)


################################################################################
###                           DATA FORMATTING                                ###
################################################################################

## Read in data file
All_CH_3yrs <- read.csv("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/NanoTag Survival/Script/CH_3yrs_final.csv")

## Sex = binary
All_CH_3yrs$Sex[All_CH_3yrs$Sex == "M"] <- 1
All_CH_3yrs$Sex[All_CH_3yrs$Sex == "F"] <- 0
All_CH_3yrs$Sex <- as.numeric(All_CH_3yrs$Sex)

## Eliminate juveniles from analysis
All_CH_3yrs<-All_CH_3yrs[!(All_CH_3yrs$Age=="FL"),]
All_CH_3yrs<-All_CH_3yrs[!(All_CH_3yrs$Age=="HY"),]

## Clean up age
All_CH_3yrs[All_CH_3yrs== "5Y"] <- "ASY"
All_CH_3yrs[All_CH_3yrs== "9Y"] <- "ASY"
All_CH_3yrs[All_CH_3yrs== "A8Y"] <- "ASY"
All_CH_3yrs[All_CH_3yrs== "AHY"] <- "ASY" 

## Age = binary
All_CH_3yrs$Age[All_CH_3yrs$Age == "SY"] <- 1
All_CH_3yrs$Age[All_CH_3yrs$Age == "ASY"] <- 0
All_CH_3yrs$Age <- as.numeric(All_CH_3yrs$Age)

## Add region column to dataset that includes all Appalachian sites (PA, TN) and
## all Great Lakes sites (WI)
All_CH_3yrs <- All_CH_3yrs %>%
  mutate(Region = paste(Site))

## Region = binary
All_CH_3yrs$Region[All_CH_3yrs$Region == "PA"] <- 0
All_CH_3yrs$Region[All_CH_3yrs$Region == "PA1"] <- 0
All_CH_3yrs$Region[All_CH_3yrs$Region == "TN"] <- 0
All_CH_3yrs$Region[All_CH_3yrs$Region == "WI"] <- 1
All_CH_3yrs$Region <- as.numeric(All_CH_3yrs$Region)

## Sites = binary
All_CH_3yrs$Site[All_CH_3yrs$Site == "PA"] <- 0
All_CH_3yrs$Site[All_CH_3yrs$Site == "PA1"] <- 1
All_CH_3yrs$Site[All_CH_3yrs$Site == "TN"] <- 2
All_CH_3yrs$Site[All_CH_3yrs$Site == "WI"] <- 3
All_CH_3yrs$Site <- as.numeric(All_CH_3yrs$Site)

## Pull in tag covariate dataframe
All_TAG <- read.csv("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/NanoTag Survival/Script/Tag_3yrs_final.csv")

## Remove non-adults from tag dataset as well
All_TAG<-All_TAG[!(All_TAG$Age=="FL"),]
All_TAG<-All_TAG[!(All_TAG$Age=="HY"),]

## Tag = binary
All_TAG[All_TAG== "1"] <- 1 #tag
All_TAG[All_TAG == "0"] <- 0 #no tag

################################################################################
###                          MODEL FORMATTING                                ###
################################################################################

## Merge CH and tag datasets together
mergedtag <- merge(All_CH_3yrs, All_TAG, by.x = "Bird.ID", by.y = "Bird.ID")
y <- as.matrix(mergedtag[,4:6])

## Get first captures
get.first <- function(x) min(which(x!=0))
f <- apply(y, 1, get.first)

## Organize covariate metadata
Tag<-mergedtag[,12]
Sex<-mergedtag[,2]
Region<-mergedtag[,8]
Site<-mergedtag[,7]
Age<-mergedtag[,3]

## Define capture history
ch <- y
colnames(ch) <- c("year1", "year2", "year3")
histories <- apply(ch[,1:3], 1, paste, collapse ="")
table(histories)

set.seed(1313)


################################################################################
###                   NIMBLE Model ONE: phi(sex), p(fixed)                   ###
################################################################################

sex.fixedp <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.sex * sex[i]
    
    gamma[1, 1, i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1, 2, i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2, 1, i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2, 2, i] <- 1           # Pr(dead t -> dead t+1)
    
  } #i
  
  
  p ~ dunif(0,1)
  
  omega[1, 1] <- 1 - p   # Pr(alive t -> non-detected t)
  omega[1, 2] <- p       # Pr(alive t -> detected t)
  omega[2, 1] <- 1           # Pr(dead t -> non-detected t)
  omega[2, 2] <- 0           # Pr(dead t -> detected t)
  
  mu ~ dnorm(mean = 0, sd = 1.5) # Prior survival intercept
  beta.sex ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter

  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu + beta.sex * sex[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2, i]) # State process
      y[i, t] ~ dcat(omega[z[i, t], 1:2]) # Observation process
    } #t
  } #i
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  p =  runif(1,0,1),
                                  z = zinits,
                                  beta.sex = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "p", "beta.sex", "mu")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
sex.fixedp <- nimbleMCMC(code = sex.fixedp,
                         constants = cjs.constants,
                         data = cjs.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = ni,
                         nburnin = nb,
                         nchains = nc,
                         thin = nt,
                         samples = TRUE,
                         WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = sex.fixedp$samples, round = 2)
sex.fixedp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = sex.fixedp$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = sex.fixedp$samples, params = c("beta.sex"))


################################################################################
###               NIMBLE Model TWO: phi(sex + tag), p(fixed)                 ###
################################################################################

sex.tag.fixedp <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.sex * sex[i] + beta.tag * tag[i]
    
    gamma[1, 1, i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1, 2, i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2, 1, i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2, 2, i] <- 1           # Pr(dead t -> dead t+1)
    
  } #i
  
  
  p ~ dunif(0,1)
  
  omega[1, 1] <- 1 - p   # Pr(alive t -> non-detected t)
  omega[1, 2] <- p       # Pr(alive t -> detected t)
  omega[2, 1] <- 1           # Pr(dead t -> non-detected t)
  omega[2, 2] <- 0           # Pr(dead t -> detected t)
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  beta.sex ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  beta.tag ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu + beta.sex * sex[i] + beta.tag * tag[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2, i]) # State process
      y[i, t] ~ dcat(omega[z[i, t], 1:2]) # Observation process
    } #t
  } #i
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex, 
                      tag = Tag)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  p =  runif(1,0,1),
                                  z = zinits,
                                  beta.sex = rnorm(1,0,1),
                                  beta.tag = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "p", "beta.sex", "beta.tag", "mu")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
sex.tag.fixedp <- nimbleMCMC(code = sex.tag.fixedp,
                             constants = cjs.constants,
                             data = cjs.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = ni,
                             nburnin = nb,
                             nchains = nc,
                             thin = nt,
                             samples = TRUE,
                             WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = sex.tag.fixedp$samples, round = 2)
sex.tag.fixedp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = sex.tag.fixedp$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = sex.tag.fixedp$samples, params = c("beta.sex", "beta.tag", 
                                                     "mean.phi", "p"))


################################################################################
###                    NIMBLE Model THREE: phi(tag), p (fixed)               ###
################################################################################

tag.fixedp <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.tag * tag[i]
    
    gamma[1, 1, i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1, 2, i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2, 1, i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2, 2, i] <- 1           # Pr(dead t -> dead t+1)
    
  } #i
  
  
  p ~ dunif(0,1)
  
  omega[1, 1] <- 1 - p   # Pr(alive t -> non-detected t)
  omega[1, 2] <- p       # Pr(alive t -> detected t)
  omega[2, 1] <- 1           # Pr(dead t -> non-detected t)
  omega[2, 2] <- 0           # Pr(dead t -> detected t)
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  beta.tag ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu + beta.tag * tag[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2, i]) # State process
      y[i, t] ~ dcat(omega[z[i, t], 1:2]) # Observation process
    } #t
  } #i
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, tag = Tag)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  p =  runif(1,0,1),
                                  z = zinits,
                                  beta.tag = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "p", "beta.tag")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
tag.fixedp <- nimbleMCMC(code = tag.fixedp,
                         constants = cjs.constants,
                         data = cjs.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = ni,
                         nburnin = nb,
                         nchains = nc,
                         thin = nt,
                         samples = TRUE,
                         WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = tag.fixedp$samples, round = 2)
tag.fixedp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = tag.fixedp$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = tag.fixedp$samples, params = c("beta.tag"))


################################################################################
###                   NIMBLE Model FOUR: phi(site), p (fixed)                ###
################################################################################

site.fixedp <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.site * site[i]
    
    gamma[1, 1, i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1, 2, i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2, 1, i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2, 2, i] <- 1           # Pr(dead t -> dead t+1)
    
  } #i
  
  
  p ~ dunif(0,1)
  
  omega[1, 1] <- 1 - p   # Pr(alive t -> non-detected t)
  omega[1, 2] <- p       # Pr(alive t -> detected t)
  omega[2, 1] <- 1           # Pr(dead t -> non-detected t)
  omega[2, 2] <- 0           # Pr(dead t -> detected t)
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  beta.site ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu + beta.site * site[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2, i]) # State process
      y[i, t] ~ dcat(omega[z[i, t], 1:2]) # Observation process
    } #t
  } #i
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, site = Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  p =  runif(1,0,1),
                                  z = zinits,
                                  beta.site = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "p", "beta.site")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
site.fixedp <- nimbleMCMC(code = site.fixedp,
                          constants = cjs.constants,
                          data = cjs.data,              
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = ni,
                          nburnin = nb,
                          nchains = nc,
                          thin = nt,
                          samples = TRUE,
                          WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = site.fixedp$samples, round = 2)
site.fixedp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = site.fixedp$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = site.fixedp$samples, params = c("beta.site"))


################################################################################
###         NIMBLE Model FIVE: phi(sex + tag + site), p (fixed)              ###
################################################################################

sex.tag.site.fixedp <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.sex * sex[i] + beta.tag * tag[i] + 
      beta.site * site[i]
    
    gamma[1, 1, i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1, 2, i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2, 1, i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2, 2, i] <- 1           # Pr(dead t -> dead t+1)
    
  } #i
  
  
  p ~ dunif(0,1)
  
  omega[1, 1] <- 1 - p   # Pr(alive t -> non-detected t)
  omega[1, 2] <- p       # Pr(alive t -> detected t)
  omega[2, 1] <- 1           # Pr(dead t -> non-detected t)
  omega[2, 2] <- 0           # Pr(dead t -> detected t)
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  beta.sex ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  beta.tag ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  beta.site ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu + beta.sex * sex[i] + beta.tag * tag[i] + 
                           beta.site * site[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2, i]) # State process
      y[i, t] ~ dcat(omega[z[i, t], 1:2]) # Observation process
    } #t
  } #i
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex, 
                      site = Site, tag = Tag)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  p =  runif(1,0,1),
                                  z = zinits,
                                  beta.sex = rnorm(1,0,1),
                                  beta.tag = rnorm(1,0,1),
                                  beta.site = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "p", "beta.sex", "beta.tag", "beta.site")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
sex.tag.site.fixedp <- nimbleMCMC(code = sex.tag.site.fixedp,
                                  constants = cjs.constants,
                                  data = cjs.data,              
                                  inits = initial.values,
                                  monitors = parameters.to.save,
                                  niter = ni,
                                  nburnin = nb,
                                  nchains = nc,
                                  thin = nt,
                                  samples = TRUE,
                                  WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = sex.tag.site.fixedp$samples, round = 2)
sex.tag.site.fixedp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = sex.tag.site.fixedp$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = sex.tag.site.fixedp$samples, params = c("beta.sex", 
                                                          "beta.tag", 
                                                          "beta.site"))


################################################################################
###                   NIMBLE Model SIX: phi(sex + site), p (fixed)           ###
################################################################################

sex.site.fixedp <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.sex * sex[i] + beta.site * site[i]
    
    gamma[1, 1, i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1, 2, i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2, 1, i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2, 2, i] <- 1           # Pr(dead t -> dead t+1)
    
  } #i
  
  
  p ~ dunif(0,1)
  
  omega[1, 1] <- 1 - p   # Pr(alive t -> non-detected t)
  omega[1, 2] <- p       # Pr(alive t -> detected t)
  omega[2, 1] <- 1           # Pr(dead t -> non-detected t)
  omega[2, 2] <- 0           # Pr(dead t -> detected t)
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  beta.sex ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  beta.site ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu + beta.sex * sex[i] + beta.site * site[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2, i]) # State process
      y[i, t] ~ dcat(omega[z[i, t], 1:2]) # Observation process
    } #t
  } #i
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex, 
                      site = Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  p =  runif(1,0,1),
                                  z = zinits,
                                  beta.sex = rnorm(1,0,1),
                                  beta.site = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "p", "beta.sex", "beta.site")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
sex.site.fixedp <- nimbleMCMC(code = sex.site.fixedp,
                              constants = cjs.constants,
                              data = cjs.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = ni,
                              nburnin = nb,
                              nchains = nc,
                              thin = nt,
                              samples = TRUE,
                              WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = sex.site.fixedp$samples, round = 2)
sex.site.fixedp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = sex.site.fixedp$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = sex.site.fixedp$samples, params = c("beta.sex",  "beta.site"))


################################################################################
###                   NIMBLE Model SEVEN: phi(fixed), p (fixed)              ###
################################################################################

intercept <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  
  
  p ~ dunif(0,1)
  phi ~ dunif(0,1)
  
  omega[1, 1] <- 1 - p   # Pr(alive t -> non-detected t)
  omega[1, 2] <- p       # Pr(alive t -> detected t)
  omega[2, 1] <- 1           # Pr(dead t -> non-detected t)
  omega[2, 2] <- 0           # Pr(dead t -> detected t)
  
  gamma[1, 1] <- phi     # Pr(alive t -> alive t+1)
  gamma[1, 2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2, 1] <- 0           # Pr(dead t -> alive t+1)
  gamma[2, 2] <- 1           # Pr(dead t -> dead t+1)
  
  
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2]) # State process
      y[i, t] ~ dcat(omega[z[i, t], 1:2]) # Observation process
    } #t
  } #i
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(p =  runif(1,0,1),
                                  phi =  runif(1,0,1),
                                  z = zinits)


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("phi", "p")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
intercept <- nimbleMCMC(code = intercept,
                        constants = cjs.constants,
                        data = cjs.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = ni,
                        nburnin = nb,
                        nchains = nc,
                        thin = nt,
                        samples = TRUE,
                        WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = intercept$samples, round = 2)
sex.site.fixedp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = intercept$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = intercept$samples, params = c("phi", "p"))


################################################################################
###                  NIMBLE Model EIGHT: phi(sex + tag), p(site)             ###
################################################################################

sex.tag.psite <- nimbleCode( {
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.sex*sex[i] + beta.tag*tag[i] #+ beta.tag*tag[i]  
    #+ beta.sex.region*region[i]*sex[i]+beta.sex.tag*sex[i]*tag[i,t]
    
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
    
    logit(p[i]) <- b0 + p.beta.site*site[i]
    
    omega[1,1,i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i]       # Pr(alive t -> detected t)
    omega[2,1,i] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0           # Pr(dead t -> detected t)
  } #i
  
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.site ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  beta.sex ~ dnorm(mean = 0, sd = 1.5)         # Prior for slope parameter
  beta.tag ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu+beta.sex*sex[i]+beta.tag*tag[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  for (i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.site*site[i])
  }
  
  mean.p <- (mean(p.est[1:2]))
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i]+1):n.occasions){
      
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:2, i]) # State process
      y[i,t] ~ dcat(omega[z[i,t], 1:2, i]) #observation process
      
    } #t
  } #i
  
  
  
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex, 
                      tag = Tag, site = Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  b0 =  rnorm(1,0,1),
                                  z = zinits, p.beta.site = rnorm(1,0,1),
                                  beta.sex = rnorm(1,0,1),
                                  beta.tag = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "mean.p", "p.beta.site", "beta.sex", 
                        "beta.tag")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
sex.tag.psite <- nimbleMCMC(code = sex.tag.psite,
                            constants = cjs.constants,
                            data = cjs.data,              
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = ni,
                            nburnin = nb,
                            nchains = nc,
                            thin = nt,
                            samples = TRUE,
                            WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = sex.tag.psite$samples, round = 2)
sex.tag.psite$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = sex.tag.psite$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)


################################################################################
###                   NIMBLE Model NINE: phi(sex), p(site)                   ###
################################################################################

sex.psite <- nimbleCode( {
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.sex*sex[i] 
    
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
    
    logit(p[i]) <- b0 + p.beta.site*site[i]
    
    omega[1,1,i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i]       # Pr(alive t -> detected t)
    omega[2,1,i] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0           # Pr(dead t -> detected t)
  } #i
  
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.site ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  beta.sex ~ dnorm(mean = 0, sd = 1.5)         # Prior for slope parameter
  
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu+beta.sex*sex[i]) 
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  for (i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.site*site[i])
  }
  
  mean.p <- (mean(p.est[1:2]))
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i]+1):n.occasions){
      
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:2, i]) # State process
      y[i,t] ~ dcat(omega[z[i,t], 1:2, i]) #observation process
      
    } #t
  } #i
  
  
  
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex, 
                      site = Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  b0 =  rnorm(1,0,1),
                                  z = zinits, p.beta.site = rnorm(1,0,1),
                                  beta.sex = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "mean.p", "p.beta.site", "beta.sex", "mu")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
sex.psite <- nimbleMCMC(code = sex.psite,
                        constants = cjs.constants,
                        data = cjs.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = ni,
                        nburnin = nb,
                        nchains = nc,
                        thin = nt,
                        samples = TRUE,
                        WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = sex.psite$samples, round = 2)
sex.psite$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = sex.psite$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)


################################################################################
###                    NIMBLE Model TEN: phi(tag), p(site)                   ###
################################################################################

tag.psite <- nimbleCode( {
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.tag*tag[i]
    
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
    
    logit(p[i]) <- b0 + p.beta.site*site[i]
    
    omega[1,1,i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i]       # Pr(alive t -> detected t)
    omega[2,1,i] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0           # Pr(dead t -> detected t)
  } #i
  
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.site ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  beta.tag ~ dnorm(mean = 0, sd = 1.5)         # Prior for slope parameter
  
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu+beta.tag*tag[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  for (i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.site*site[i])
    
  }
  
  mean.p <- (mean(p.est[1:2]))
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i]+1):n.occasions){
      
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:2, i]) # State process
      y[i,t] ~ dcat(omega[z[i,t], 1:2, i]) #observation process
      
    } #t
  } #i
  
  
  
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, tag = Tag, 
                      site = Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  b0 =  rnorm(1,0,1),
                                  z = zinits, p.beta.site = rnorm(1,0,1),
                                  beta.tag = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "mean.p", "p.beta.site", "beta.tag")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
tag.psite <- nimbleMCMC(code = tag.psite,
                        constants = cjs.constants,
                        data = cjs.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = ni,
                        nburnin = nb,
                        nchains = nc,
                        thin = nt,
                        samples = TRUE,
                        WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = tag.psite$samples, round = 2)
tag.psite$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = tag.psite$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)


################################################################################
###                 NIMBLE Model ELEVEN: phi(site), p(site)                  ###
################################################################################

site.psite <- nimbleCode( {
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.site*site[i] 
    
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
    
    logit(p[i]) <- b0 + p.beta.site*site[i]
    
    omega[1,1,i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i]       # Pr(alive t -> detected t)
    omega[2,1,i] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0           # Pr(dead t -> detected t)
  } #i
  
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.site ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  beta.site ~ dnorm(mean = 0, sd = 1.5)         # Prior for slope parameter
  
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu+beta.site*site[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  for (i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.site*site[i])
    
  }
  
  mean.p <- (mean(p.est[1:2]))
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i]+1):n.occasions){
      
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:2, i]) # State process
      y[i,t] ~ dcat(omega[z[i,t], 1:2, i]) #observation process
      
    } #t
  } #i
  
  
  
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, site = Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  b0 =  rnorm(1,0,1),
                                  z = zinits, p.beta.site = rnorm(1,0,1),
                                  beta.site = rnorm(1,0,1))

######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mean.phi", "mean.p", "p.beta.site", "beta.site")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
site.psite <- nimbleMCMC(code = site.psite,
                         constants = cjs.constants,
                         data = cjs.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = ni,
                         nburnin = nb,
                         nchains = nc,
                         thin = nt,
                         samples = TRUE,
                         WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = site.psite$samples, round = 2)
site.psite$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = site.psite$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)


################################################################################
###             NIMBLE Model TWELVE: phi(sex + site + tag), p(site)          ###
################################################################################

sex.site.tag.psite <- nimbleCode( {
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.sex*sex[i] + beta.site*site[i] + beta.tag*tag[i]
    
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
    
    logit(p[i]) <- b0 + p.beta.site*site[i]
    
    omega[1,1,i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i]       # Pr(alive t -> detected t)
    omega[2,1,i] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0           # Pr(dead t -> detected t)
  } #i
  
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.site ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  beta.sex ~ dnorm(mean = 0, sd = 1.5)         # Prior for slope parameter
  beta.site ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  beta.tag ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu+beta.sex*sex[i]+beta.site*site[i]+beta.tag*tag[i])
  }
  
  mean.phi <- (mean(phi.est[1:nind]))
  
  
  for (i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.site*site[i])
    
  }
  
  mean.p <- (mean(p.est[1:nind]))
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i]+1):n.occasions){
      
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:2, i]) # State process
      y[i,t] ~ dcat(omega[z[i,t], 1:2, i]) #observation process
      
    } #t
  } #i
  
  
  
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex, 
                      site = Site, tag = Tag)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  b0 =  rnorm(1,0,1),
                                  z = zinits, 
                                  p.beta.site = rnorm(1,0,1),
                                  beta.sex = rnorm(1,0,1),
                                  beta.site = rnorm(1,0,1),
                                  beta.tag = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("phi.est", "mean.phi", "mean.p", "p.beta.site", 
                        "beta.sex", "beta.site", "beta.tag", "mu", "b0")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
sex.site.tag.psite <- nimbleMCMC(code = sex.site.tag.psite,
                                 constants = cjs.constants,
                                 data = cjs.data,              
                                 inits = initial.values,
                                 monitors = parameters.to.save,
                                 niter = ni,
                                 nburnin = nb,
                                 nchains = nc,
                                 thin = nt,
                                 samples = TRUE,
                                 WAIC = TRUE)

#######pull out estimate summary, check Rhat
sex.site.tag.psite(object = sex.site.tag.psite$samples, round = 2)
sex.site.tag.psite$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = sex.site.tag.psite$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)


beta.site = rnorm(1,0,1)

MCMCsummary(object = sex.site.tag.psite$samples, round = 2)
MCMCplot(object = sex.site.tag.psite$samples, params = c("beta.sex",  "beta.site", "beta.tag"))


################################################################################
###             NIMBLE Model THIRTEEN: phi(sex + site), p(site)              ###
################################################################################

sex.site.psite <- nimbleCode( {
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  for (i in 1:nind){
    
    logit(phi[i]) <- mu + beta.sex*sex[i] + beta.site*site[i]
    
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
    
    logit(p[i]) <- b0 + p.beta.site*site[i]
    
    omega[1,1,i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i]       # Pr(alive t -> detected t)
    omega[2,1,i] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0           # Pr(dead t -> detected t)
  } #i
  
  
  mu ~ dnorm(mean = 0, sd = 1.5) #prior survival intercept
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.site ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  beta.sex ~ dnorm(mean = 0, sd = 1.5)         # Prior for slope parameter
  beta.site ~ dnorm(mean = 0, sd = 1.5) # Prior for slope parameter
  
  for (i in 1:nind){
    phi.est[i] <- ilogit(mu+beta.sex*sex[i]+beta.site*site[i])
  }
  
  mean.phi <- (mean(phi.est[1:2]))
  
  
  for (i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.site*site[i])
    
  }
  
  mean.p <- (mean(p.est[1:2]))
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i]+1):n.occasions){
      
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:2, i]) # State process
      y[i,t] ~ dcat(omega[z[i,t], 1:2, i]) #observation process
      
    } #t
  } #i
  
  
  
})

cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex, 
                        site = Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(mu =  rnorm(1,0,1),
                                  b0 =  rnorm(1,0,1),
                                  z = zinits, p.beta.site = rnorm(1,0,1),
                                  beta.sex = rnorm(1,0,1),
                                  beta.site = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("mu", "b0", "mean.phi", "mean.p", "p.beta.site", 
                        "beta.sex", "beta.site")
## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
sex.site.psite <- nimbleMCMC(code = sex.site.psite,
                             constants = cjs.constants,
                             data = cjs.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = ni,
                             nburnin = nb,
                             nchains = nc,
                             thin = nt,
                             samples = TRUE,
                             WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = sex.site.psite$samples, round = 2)
sex.site.psite$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = sex.site.psite$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)


################################################################################
###                 NIMBLE Model FOURTEEN: phi(fixed), p (site)              ###
################################################################################

intercept.sitep <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  
  # Priors and constraints
  for (i in 1:nind){
    
    logit(p[i]) <- b0 + p.beta.site*site[i]
    
    omega[1,1,i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i]       # Pr(alive t -> detected t)
    omega[2,1,i] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0           # Pr(dead t -> detected t)
  } #i
  
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.site ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  
  phi ~ dunif(0,1)
  
  gamma[1, 1] <- phi     # Pr(alive t -> alive t+1)
  gamma[1, 2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2, 1] <- 0           # Pr(dead t -> alive t+1)
  gamma[2, 2] <- 1           # Pr(dead t -> dead t+1)
  
  
  for (i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.site*site[i])
    
  }
  
  mean.p <- (mean(p.est[1:2]))  
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] ~ dcat(delta[1:2])
    
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2]) # State process
      y[i, t] ~ dcat(omega[z[i, t], 1:2, i]) # Observation process
    } #t
  } #i
})



cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, site=Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi =  runif(1,0,1),
                                  z = zinits,
                                  b0 =  rnorm(1,0,1),
                                  p.beta.site = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("phi", "mean.p", "p.beta.site")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
intercept.sitep <- nimbleMCMC(code = intercept.sitep,
                              constants = cjs.constants,
                              data = cjs.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = ni,
                              nburnin = nb,
                              nchains = nc,
                              thin = nt,
                              samples = TRUE,
                              WAIC = TRUE)

#######pull out estimate summary, check Rhat
MCMCsummary(object = intercept$samples, round = 2)
intercept$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = intercept$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = intercept$samples, params = c("phi", "p"))


################################################################################
###                       WAIC TO DETERMINE TOP MODEL                        ###
################################################################################

## Create WAIC table
WAIC_scores <- data.frame(model = c("(sex.tag.psite)",
                                    "(sex.psite)",
                                    "(tag.psite)",
                                    "(site.psite)",
                                    "(sex.site.tag.psite)",
                                    "(sex.site.psite)",
                                    "(intercept)",
                                    "(sex.tag.fixedp)",
                                    "(sex.fixedp)",
                                    "(tag.fixedp)",
                                    "(site.fixedp)",
                                    "(sex.tag.site.fixedp)",
                                    "(sex.site.fixedp)",
                                    "(intercept.sitep)"),
                          WAIC = c(sex.tag.psite$WAIC$WAIC,
                                   sex.psite$WAIC$WAIC,
                                   tag.psite$WAIC$WAIC,
                                   site.psite$WAIC$WAIC,
                                   sex.site.tag.psite$WAIC$WAIC,
                                   sex.site.psite$WAIC$WAIC,
                                   intercept$WAIC$WAIC,
                                   sex.tag.fixedp$WAIC$WAIC,
                                   sex.fixedp$WAIC$WAIC,
                                   tag.fixedp$WAIC$WAIC,
                                   site.fixedp$WAIC$WAIC,
                                   sex.tag.site.fixedp$WAIC$WAIC,
                                   sex.site.fixedp$WAIC$WAIC,
                                   intercept.sitep$WAIC$WAIC))

## Plot beta coefficients
MCMCplot(object = tag.fixedp$samples, params = c("beta.tag"))
MCMCplot(object = sex.fixedp$samples, params = c("beta.sex"))
MCMCplot(object = site.fixedp$samples, params = c("beta.site"))
MCMCplot(object = region.fixedp$samples, params = c("beta.region"))
MCMCplot(object = sex.tag.site.fixedp$samples, params = c("beta.sex",
                                                          "beta.site",
                                                          "beta.tag"))
