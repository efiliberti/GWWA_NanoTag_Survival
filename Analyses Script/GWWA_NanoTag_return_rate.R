## This script runs a series of GLMMs to identify relationships between return
## rate (predictor variable) and spatial and demographic covariates (response
## variables) of Golden-winged Warbler. CMR datasets collected between 2021-2023
## and include both 3-year datasets and 2-year datasets. 3-year individuals that 
## were reencountered in the second occastion were included as two different 
## records (*technically* different birds) representing 2021-2022 and 2022-2023.

## Load packages
library(ggplot2)
library(lme4)
library(sjPlot)
library(dplyr)
library(performance)
library(cowplot)
library(ggeffects)

################################################################################
###                       FORMAT CAPTURE HISTORY                             ###
################################################################################

## Read in NanoTag capture history
CH <- read.csv("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/NanoTag Survival/Script/CH_2yrs_final.csv")

## Put search effort (hrs per site) into dataframe
search <- data.frame(Effort = as.integer(c(95.35, 190.75, 34.75, 34.75, 65.90, 
                                           14.75, 45.05, 18.86, 10.10, 17.20, 
                                           17.20, 4.00, 629.07, 294.50, 34.75, 
                                           34.75)),
                     Site = c("TN", "WI", "PA", "PA1", "VA", "MNB", "NCK", "NY", 
                              "KY", "NC", "NC1", "MN", "TN", "WI", "PA", "PA1"),
                     Year.Banded = c(2022, 2022, 2022, 2022, 2022, 2022, 2022, 
                              2022, 2022, 2022, 2022, 2022, 2021, 2021, 2021, 
                              2021))

## Merge search effort with capture history
CH <- left_join(CH, search, by = c("Site", "Year.Banded"))

## Standardize effort data 
CH$Effort <- datawizard::standardize((CH$Effort))

## Remove hatch-year individuals from the analysis
CH<-CH[!(CH$Age=="FL"),]
CH<-CH[!(CH$Age=="HY"),]

## Change all after-hatch-year (AHY) individuals to after-second-year (ASY)
## individuals (n = 13) so that you can run the most parsimonious model without 
## sacrificing sample size for non-age models. 
CH[CH== "AHY"] <- "ASY" 

## Convert individual sites into respective regions.
CH$Region[CH$Site == "PA"] <- "App"
CH$Region[CH$Site == "PA1"] <- "App"
CH$Region[CH$Site == "TN"] <- "App"
CH$Region[CH$Site == "NC"] <- "App"
CH$Region[CH$Site == "NC1"] <- "App"
CH$Region[CH$Site == "NCK"] <- "App"
CH$Region[CH$Site == "VA"] <- "App"
CH$Region[CH$Site == "KY"] <- "App"
CH$Region[CH$Site == "WI"] <- "GL"
CH$Region[CH$Site == "MN"] <- "GL"
CH$Region[CH$Site == "MNB"] <- "GL"
CH$Region[CH$Site == "NY"] <- "GL"

## Rename sites for easier analysis
CH$Site[CH$Site == "PA"] <- "Pennsylvania (E)"
CH$Site[CH$Site == "PA1"] <- "Pennsylvania (C)"
CH$Site[CH$Site == "TN"] <- "Tennessee"
CH$Site[CH$Site == "NC1"] <- "North Carolina (NW)"
CH$Site[CH$Site == "NC"] <- "North Carolina (W)"
CH$Site[CH$Site == "NCK"] <- "North Carolina (SW)"
CH$Site[CH$Site == "VA"] <- "Virginia"
CH$Site[CH$Site == "KY"] <- "Kentucky"
CH$Site[CH$Site == "WI"] <- "Wisconsin"
CH$Site[CH$Site == "MN"] <- "Minnesota (C)"
CH$Site[CH$Site == "MNB"] <- "Minnesota (N)"
CH$Site[CH$Site == "NY"] <- "New York"

## Order site levels and convert to factor
order <- c("Wisconsin", "New York", "Minnesota (C)", "Minnesota (N)", 
           "Tennessee", "Pennsylvania (E)", "Pennsylvania (C)", "Kentucky", 
           "North Carolina (SW)", "North Carolina (W)", "North Carolina (NW)", 
           "Virginia")
all_CH$Site <- factor(all_CH$Site, levels = order)

## Rename cohort for easier analysis
CH$Tag[CH$Cohort == "T"] <- "Present"
CH$Tag[CH$Cohort == "C"] <- "Absent"

## Year banded = categorical
CH$Year.Banded[CH$Year.Banded == 2021] <- "2021"
CH$Year.Banded[CH$Year.Banded == 2022] <- "2022"

## Collect metadata for dataset
all_CH <- CH
all_CH$sum <- rowSums(CH[,6:7])
table(all_CH$CH)
table(all_CH$sum)
all_CH$Return <- ifelse(all_CH$sum >1, 1, 0)
table(all_CH$Return)

set.seed(1313)

################################################################################
###                                   GLMMs                                  ###
################################################################################

library(lme4)
library(AICcmodavg)

## GLMM with site included as a random effect -- various models looking at 
## relationships among sex, cohort, region, year banded, age, and hours spent 
## searching per site.

m.Intercept <- glmer(Return ~ 1 + (1|Site), data = all_CH, family = binomial)
m.Sex <- glmer(Return ~ Sex + (1|Site), data = all_CH, family = binomial)
m.Tag <- glmer(Return ~ Tag + (1|Site), data = all_CH, family = binomial)
m.Region <- glmer(Return ~ Region + (1|Site), data = all_CH, family = binomial)
m.Year <- glmer(Return ~ Year.Banded + (1|Site), data = all_CH, family = 
                  binomial)
m.Sex.a.Tag <- glmer(Return ~ Sex + Tag + (1|Site), data = all_CH, family = 
                       binomial)
m.Age <- glmer(Return ~ Age + (1|Site), data = all_CH, family = binomial)
m.Sex.t.Region <- glmer(Return ~ Sex * Region + (1|Site), data = all_CH, family 
                        = binomial)
m.Effort <- glmer(Return ~ Effort + (1|Site), data = all_CH, family = binomial)
m.StudyCovs <- glmer(Return ~ Effort + Year.Banded + Tag + (1|Site), data = 
                       all_CH, family = binomial)
m.DemCovs <- glmer(Return ~ Sex + Age + Region + (1|Site), data = all_CH, 
                   family = binomial)


## Pull into an AIC table for model comparison
AIC <- aictab(list(m.Intercept, m.Sex, m.Tag, m.Region, m.Year, m.Age, 
                   m.Sex.t.Region, m.Effort, m.StudyCovs, m.DemCovs, 
                   m.Sex.a.Tag), modnames = c("Intercept", "Sex", "Tag",  
                                              "Region", "Year", "Age", 
                                              "Sex*Region",  "Effort", 
                                              "all.study.cov", "all.cov", 
                                              "Sex.Tag"),
              dname = "all_CH")


# Find the mean return rate for all individuals (and SD) for top performing 
## model (sex as a fixed effect and site as a random effect).
rr.est <- predict(m.Sex, type = "response")
rr.phi <- mean(rr.est)
rr.sd <- sd(rr.est)
rr.low.ci <- quantile(rr.est, probs = c(2.5)/100)
rr.upp.ci <- quantile(rr.est, probs = c(97.5)/100)

## Obtain female return rate
female_data <- all_CH %>%
  filter(Sex == "F")
rr_est_f <- predict(m.Sex, newdata = female_data, type = "response")
rr_est_f_mean <- mean(rr_est_f)
rr_est_f_sd <- sd(rr_est_f)
rr_up_f <- mean(rr_est_f) + sd(rr_est_f)
rr_lo_f <- mean(rr_est_f) - sd(rr_est_f)

## Obtain male return rate
male_data <- all_CH %>%
  filter(Sex == "M")
rr_est_m <- predict(m.Sex, newdata = male_data, type = "response")
rr_est_m_mean <- mean(rr_est_m)
rr_est_m_sd <- sd(rr_est_m)
rr_up_m <- mean(rr_est_m) + sd(rr_est_m)
rr_lo_m <- mean(rr_est_m) - sd(rr_est_m)
