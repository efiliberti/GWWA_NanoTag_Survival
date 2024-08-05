
## Figure 3 in manuscript. This visualization pulls estimates out of two models: 
## one that includes tag, year, and effort (study design covariates) as additive 
## fixed effects and one that includes sex, region, and age (spatial/demographic 
## covariates) as additive fixed effects. Site is included as a random effect 
## for all models.

## Load packages
library(ggplot2)
library(lme4)
library(sjPlot)
library(dplyr)
library(performance)
library(cowplot)


################################################################################
####                         ODDS RATIOS OF COVARIATES                      ####
################################################################################

plot1 <- plot_model(m.StudyCovs, type = "est", show.values = T, 
                    title = "", colors = "black", grid = T, axis.lim = c(0.3, 2.9))+
  theme_bw() + 
  geom_hline(yintercept=1, linetype="dashed", color = "red", linewidth=0.5) +
  ggtitle("Study Design Covariates")+
  ggeasy::easy_center_title()+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "transparent"),
        panel.grid.minor = element_line(color = "transparent"),
        axis.text.x = element_text(color = "black", size = 14), 
        axis.text.y = element_text(color = "black", size = 14),  
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 14.5))

setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Thesis/Chapter 3/Figures/Figure 3.3")
ggsave("plot1.png", dpi = 1200, width = 4, height = 5, units = "in")


plot2 <- plot_model(m.DemCovs, type = "est", show.values = T, 
                    title = "", colors = "black", grid = T, axis.lim = c(0.3, 2.9))+
  theme_bw() + 
  geom_hline(yintercept=1, linetype="dashed", color = "red", linewidth=0.5) +
  ggtitle("Spatial/Demographic Covariates")+
  ggeasy::easy_center_title()+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "transparent"),
        panel.grid.minor = element_line(color = "transparent"),
        axis.text.x = element_text(color = "black", size = 14),  
        axis.text.y = element_text(color = "black", size = 14), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 14.5))

ggsave("plot2.png", dpi = 1200, width = 4, height = 5, units = "in")


plot_grid(plot1, plot2)
