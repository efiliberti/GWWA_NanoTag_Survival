
## Figure A4 in manuscript. This script shows site-specific random effect 
## contribution to explore variations in return rates among sites. Pulled from 
## sex as a solo fixed effect. 

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
####                            RANDOM SITE EFFECT                          ####
################################################################################

## Random effect of site for selected model
site.random.effects <- plot_model(m.Sex, type = "re", color = "black", grid = T, show.p = T)+
  theme_bw()+
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
        legend.position = "none")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  geom_vline(xintercept = 4.5, size = 1, color = "gray")+
  ylim(0, 2.5)+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Odds Ratio")


setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/NanoTag Survival/Tables and Figures/Figure A4")
ggsave("site.random.effects.png", dpi = 1200, width = 4, height = 5, units = "in")

