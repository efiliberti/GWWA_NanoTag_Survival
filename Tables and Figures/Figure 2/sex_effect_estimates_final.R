
## Figure 2 in manuscript. This visualization pulls and compares sex effect 
## estimates from the annual survival and return rate analyses. Fixed effect of 
## sex for return rate model and fixed effects of sex, site, and tag for annual 
## survival analysis. 

## Load packages
library(ggplot2)
library(lme4)
library(sjPlot)
library(dplyr)
library(performance)
library(cowplot)

########################## ANNUAL SURVIVAL ESTIMATE ##########################

# Extract samples for mean.phi.f and mean.phi.m
as_est_f <- sex.site.tag.psite$samples$chain1[, 'mean.phi.f']
as_est_f <- c(as_est_f, sex.site.tag.psite$samples$chain2[, 'mean.phi.f'])

as_est_m <- sex.site.tag.psite$samples$chain1[, 'mean.phi.m']
as_est_m <- c(as_est_m, sex.site.tag.psite$samples$chain2[, 'mean.phi.m'])

# Calculate credible intervals
as_up_f <- quantile(as_est_f, probs = 0.975)
as_lo_f <- quantile(as_est_f, probs = 0.025)
as_up_m <- quantile(as_est_f, probs = 0.975)
as_lo_m <- quantile(as_est_f, probs = 0.025)

## Put metadata into dataframe
as_sex <- data.frame(Sex = c("Female", "Male"),
                     Mean = c(as_est_f, as_est_m),
                     Lower = c(as_lo_f, as_lo_m),
                     Upper = c(as_up_f, as_up_m),
                     Analysis = c("Annual Survival", "Annual Survival"))


############################ RETURN RATE ESTIMATES #############################

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

## Put metadata into dataframe
rr_sex <- data.frame(Sex = c("Female", "Male"),
                     Mean = c(rr_est_f_mean, rr_est_m_mean),
                     Lower = c(rr_lo_f, rr_lo_m),
                     Upper = c(rr_up_f, rr_up_m),
                     Analysis = c("Annual Return Rate", "Annual Return Rate"))


############################### COMBINE DATAFRAMES #############################

as_rr_viz <- rbind(as_sex, rr_sex)
as_rr_viz$Analysis <- factor(as_rr_viz$Analysis, levels = 
                               c("Annual Return Rate", 
                                 "Apparent Annual Survival"))


############################### DATA VISUALIZATION #############################

## Create plot
as_rr_viz_plot <- ggplot(as_rr_viz, aes(y = Mean, x = factor(Sex), shape = factor(Analysis, levels = rev(levels(Analysis))), fill = factor(Analysis, levels = rev(levels(Analysis))))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.05, position = position_dodge(width = 0.4)) +
  geom_point(size = 2, position = position_dodge(width = 0.4), color = "black") +
  theme(panel.grid.major = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 13),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title = element_text(colour = "white"),
        legend.position = "none") +
  ylab("Annual Survival and Return Rate") +
  xlab("") +
  theme(panel.grid.minor = element_blank()) +
  ylim(0, 0.5) +
  scale_shape_manual(values = c(21, 24)) +  # 21 for circles, 24 for triangles
  scale_fill_manual(values = c("black", "white"))+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "transparent"),
        panel.grid.minor = element_line(color = "transparent")
  )


setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/NanoTag Survival/Tables and Figures/Figure 2")
ggsave("as_rr_viz_plot.png", dpi = 1200, width = 2, height = 3, units = "in")
