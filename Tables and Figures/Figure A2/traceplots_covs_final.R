
## Figure A2 in manuscript. This script runs the selected phi(sex+site+tag), 
## p(site) apparent annual survival model that was selected to interpret beta 
## coefficients from. Traceplots will show siccessfulconvergence of model. 

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

################################################################################
###                                TRACEPLOTS                                ###
################################################################################

setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/NanoTag Survival/Tables and Figures/Figure A2/")

traceplot <- MCMCtrace(object = sex.site.tag.psite$samples,
          ISB = FALSE,
          exact = TRUE,
          pdf = TRUE)

dev.off()
