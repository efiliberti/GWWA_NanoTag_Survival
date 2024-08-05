
## Figure A1 in manuscript. This visualization shows standardized total 
## search effort (hr) spent looking for control and NanoTag birds. Raw total 
## search effort hours were pulled from datasheets contributed by each 
## collaborator. 

## Load packages
library(ggplot2)
library(lme4)
library(sjPlot)
library(dplyr)
library(performance)
library(cowplot)


################################################################################
####                     LOG-TRANSFORMED SERACH EFFORT                      ####
################################################################################

search_effort <- data.frame(Site = c("Wisconsin", "Wisconsin", "New York", "Minnesota (C)",
                                      "Minnesota (N)", "Tennessee", "Tennessee",
                                      "Pennsylvania", "Pennsylvania", "Kentucky", "North Carolina (N)",
                                      "North Carolina (SW)", "Virginia"),
                             Year = c("2021", "2022","2022", "2022", "2022", "2022", "2021", "2021", "2022",
                                      "2022","2022","2022", "2022"),
                            Hour = c(294.50, 190.75, 18.86, 14.75, 4, 95.35, 629.07, 34.75, 34.75, 10.10, 
                                     45.05, 17.20, 65.90))

search_effort$Site <- factor(search_effort$Site, levels = unique(search_effort$Site))
search_effort$logHour <- datawizard::standardize((search_effort$Hour))

ggplot(data = search_effort, aes(x = Site, y = logHour, fill = Year)) +
  geom_col(width = 0.75)+
  geom_vline(xintercept=4.5, linetype="dashed", color = "red", size=0.5) +
  theme_bw()+
  xlab("Site")+ ylab("logHour")+
  scale_fill_manual(values = c('darkgray', 'black'))+
  theme(axis.text.x = element_text(angle = 60, hjust = 01))+
  theme(panel.grid.minor = element_blank())+
  ggtitle("NanoTag Search Effort")+
  theme(plot.title = element_text(hjust = 0.5))
