## This script runs apparent annual survival CJS models on Golden-winged Warbler 
## CMR datasets collected between 2021-2023. We test different spatial and 
## demographic covariates as fixed effects on detection probability to see what 
## model explains the most variation. This is the first modelset we created for
## our manuscript survival analysis. 

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
###                   NIMBLE Model ONE: phi(fixed), p (tag)                  ###
################################################################################


## Model that includes sex and region as fixed effects on survival and a tag
## effect on detection probability

fixedphi.tagp <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  
  phi ~ dunif(0,1)
  
  gamma[1, 1] <- phi     # Pr(alive t -> alive t+1)
  gamma[1, 2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2, 1] <- 0           # Pr(dead t -> alive t+1)
  gamma[2, 2] <- 1           # Pr(dead t -> dead t+1)
  
  for (i in 1:nind){
    
    logit(p[i]) <- b0 + p.beta.tag*tag[i]
    
  omega[1, 1, i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
  omega[1, 2, i] <- p[i]       # Pr(alive t -> detected t)
  omega[2, 1, i] <- 1           # Pr(dead t -> non-detected t)
  omega[2, 2, i] <- 0           # Pr(dead t -> detected t)
  
  } #i
  
  
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.tag ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  
  for(i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.tag*tag[i])
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


cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, tag = Tag)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(b0 =  runif(1,0,1),
                                  phi = runif(1,0,1),
                                  z = zinits,
                                  p.beta.tag = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("phi", "mean.p", "p.beta.tag")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
fixedphi.tagp <- nimbleMCMC(code = fixedphi.tagp,
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
MCMCsummary(object = fixedphi.tagp$samples, round = 2)
fixedphi.tagp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = fixedphi.tagp$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = fixedphi.tagp$samples, params = c("p.beta.tag"))


################################################################################
###                   NIMBLE Model TWO: phi(fixed), p (sex)                  ###
################################################################################

fixedphi.sexp <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  
  phi ~ dunif(0,1)
  
  gamma[1, 1] <- phi     # Pr(alive t -> alive t+1)
  gamma[1, 2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2, 1] <- 0           # Pr(dead t -> alive t+1)
  gamma[2, 2] <- 1           # Pr(dead t -> dead t+1)
  
  for (i in 1:nind){
    
    logit(p[i]) <- b0 + p.beta.sex*sex[i]
    
    omega[1, 1, i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1, 2, i] <- p[i]       # Pr(alive t -> detected t)
    omega[2, 1, i] <- 1           # Pr(dead t -> non-detected t)
    omega[2, 2, i] <- 0           # Pr(dead t -> detected t)
    
  } #i
  
  
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.sex ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  
  for(i in 1:nind){
    p.est[i] <- ilogit(b0 + p.beta.sex*sex[i])
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


cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, sex = Sex)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(b0 =  runif(1,0,1),
                                  phi = runif(1,0,1),
                                  z = zinits,
                                  p.beta.sex = rnorm(1,0,1))


######### Bundle data

## Parameters monitored ########
parameters.to.save <- c("phi", "mean.p", "p.beta.sex")

## MCMC settings
ni <- 50000
nt <- 3
nb <- 10000
nc <- 3

##################### Run and summarize the model ##########################
fixedphi.sexp <- nimbleMCMC(code = fixedphi.sexp,
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
MCMCsummary(object = fixedphi.sexp$samples, round = 2)
fixedphi.sexp$WAIC

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = fixedphi.sexp$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = fixedphi.sexp$samples, params = c("phi", "mean.p", 
                                                    "p.beta.sex"))


################################################################################
###                   NIMBLE Model THREE: phi(fixed), p (site)               ###
################################################################################

fixedphi.sitep <- nimbleCode({
  #######################################################################
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # Priors and constraints
  
  phi ~ dunif(0,1)
  
  gamma[1, 1] <- phi     # Pr(alive t -> alive t+1)
  gamma[1, 2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2, 1] <- 0           # Pr(dead t -> alive t+1)
  gamma[2, 2] <- 1           # Pr(dead t -> dead t+1)
  
  for (i in 1:nind){
    
    logit(p[i]) <- b0 + p.beta.site*site[i]
    
    omega[1, 1, i] <- 1 - p[i]   # Pr(alive t -> non-detected t)
    omega[1, 2, i] <- p[i]       # Pr(alive t -> detected t)
    omega[2, 1, i] <- 1           # Pr(dead t -> non-detected t)
    omega[2, 2, i] <- 0           # Pr(dead t -> detected t)
    
  } #i
  
  
  b0 ~ dnorm(mean = 0, sd = 1.5) #prior recapture intercept
  p.beta.site ~ dnorm(mean = 0, sd = 1.5)        # Prior for slope parameter
  
  for(i in 1:nind){
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


cjs.constants <- list(nind=nrow(y), n.occasions=ncol(y), f=f, site = Site)
cjs.data <- list(y=(y+1))


zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(b0 =  runif(1,0,1),
                                  phi = runif(1,0,1),
                                  z = zinits,
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
fixedphi.sitep <- nimbleMCMC(code = fixedphi.sitep,
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
MCMCsummary(object = fixedphi.sitep$samples, round = 2)


###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = fixedphi.sitep$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = fixedphi.sitep$samples, params = c("phi", "mean.p", "p.beta.site"))


################################################################################
###                    NIMBLE Model FOUR: phi(fixed), p (fixed)              ###
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

###############Trace plots to check parameter mixing and chain convergence
MCMCtrace(object = intercept$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCplot(object = intercept$samples, params = c("phi", "p"))



################################################################################
###                        WAIC TO DETERMINE TOP MODEL                       ###
################################################################################

## Create WAIC table
p_WAIC_scores <- data.frame(model = c("(fixedphi.tagp)",
                                    "(fixedphi.sexp)",
                                    "(fixedphi.sitep)",
                                    "(intercept)"),
                          WAIC = c(fixedphi.tagp$WAIC$WAIC,
                                   fixedphi.sexp$WAIC$WAIC,
                                   fixedphi.sitep$WAIC$WAIC,
                                   intercept$WAIC$WAIC))
