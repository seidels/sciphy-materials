## ---------------------------
##
## Script name: plot_growth_rates_OU
##
##
## Author: Sophie Seidel
##
## Date Created: 2024-05-13
##
## Copyright (c) Sophie Seidel 2024
## Email: sophie.seidel@posteo.de

## set working directory where the log files are
setwd("~/Projects/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target_OU/")
figure_dir = "~/Projects/sciphy-materials/figure_2/"

# text size
text_size = 10

## load up the packages we will need:
library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)
library(pammtools)
library(reshape2)


timeline_format <- c(seq(0,25,by=2), 25)

# Define helper functions
get_median_and_hpd = function(growth){

  HPD <- HPDinterval(growth)

  name = paste0("growthRate.", 1:ncol(growth))
  median <- as.numeric(sapply( data.frame(growth),median))
  up_bd <- as.numeric(HPD[,"upper"])
  low_bd <- as.numeric(HPD[,"lower"])

  return(data.frame(name=name, median=median, hpd_up = up_bd, hpd_low=low_bd))
}
get_determ_pop_size = function(median){
  total <- 1
  for(i in 1:12) {

    total <- total*exp(median[i]*2)

  }
  total <- total*exp(median[13])

  return(total)
}

#load the log file
typewriter_file <- "combined_burnin10.log"
typewriter <- read.table(typewriter_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] = c("growthRate.14", growth_hpd[13, 2:4]) # add last timepoint, st for time 24 - 25 the same value remains
get_determ_pop_size(median = growth_hpd$median) # 848055.3

growth_hpd$t = timeline_format

p_growth = ggplot(growth_hpd, aes(x=t, y=median)) +
  #facet_grid(rows =  (startsWith(x = tree, prefix = "Sci")))+
  #facet_grid(tree~.)+
  #geom_line() +
  geom_step(size=1)+
  geom_stepribbon(aes(ymin = hpd_low, ymax=hpd_up), linetype="dotted", alpha = 0.1)+
  #geom_ribbon(aes(ymin = hpd_low, ymax=hpd_up, fill=tree), linetype="dotted", alpha = 0.1)+
  theme_classic() +
  theme(legend.title=element_blank(), legend.position = c(0.27, 0.858), text = element_text(size = text_size)) +
  ylab(expression("Dynamic growth rate [" * d^-1 * "]"))+
  xlab("Time [d]")
p_growth

ggsave(paste0(figure_dir, "growth_skyline_ou.png"), p_growth , width = 9.52, height = 4.76, units = "cm", dpi = 300)
svg(filename = paste0(figure_dir, "growth_skyline_ou.svg"), width = 9.52, height = 4.76, pointsize = text_size)
p_growth
dev.off()

