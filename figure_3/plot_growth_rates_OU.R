## ---------------------------
##
## Script name: plot_growth_rates
##
## Purpose of script: Plot estimates from SciPhy on HEK293 cell culture data, skyline clock per target, with color annotations matching the tree.
##
## Author: Antoine Zwaans
##
## Date Created: 2023-10-13
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch

## set working directory where the log files are
setwd("~/Projects/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target_OU/")
figure_dir = "~/Projects/sciphy-materials/figure_2/"

## load up the packages we will need:
library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)

#load the log file
typewriter_file <- "combined_burnin10.log"
typewriter <- read.table(typewriter_file, header = T)

# text size
text_size = 10
#########################################
##plotting growth rates as violin plots##
#########################################

birth_rates <- typewriter[,paste0(rep("birthRate.",13),1:13)]
death_rates <- typewriter[,paste0(rep("deathRate.",13),1:13)]

growth_rates <- birth_rates - death_rates
##add a column with prior
growth_rates <- bind_cols(growth_rates,Prior=rlnorm(nrow(growth_rates), meanlog = 0.1, sdlog = 1.0) - rlnorm(nrow(growth_rates), meanlog = -0.4, sdlog = 1.0))
growth_rates <- growth_rates[order(-sapply(growth_rates, median))]



names(growth_rates) <- c(paste0(rep("growthRate.",13),1:13),"Prior")
growth_rates_long <- pivot_longer(growth_rates,seq(1,ncol(growth_rates)))
ordered_names <- c(paste0(rep("growthRate.",13),1:13))
growth_rates_long <- mutate(growth_rates_long,name = fct_relevel(name,ordered_names,"Prior"))

p <- ggplot(growth_rates_long,aes(x=name,value,fill=name)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Growth") +
  ylab(parse(text = paste0('"Posterior growth rate "', '(~d^-1)'))) +
  scale_fill_manual(values=c(rep("#5CA17D",13),"#E1E1F7")) +
  coord_cartesian(ylim = c(-1.0,4.0),expand = TRUE) +
  theme(legend.position = "none" ,axis.text.x = element_text(angle = 90)) +
  theme(text=element_text(size = text_size),panel.grid.minor = element_blank(),
                                                                                    panel.border = element_blank(),
                                                                                    panel.background = element_blank(),panel.grid.major.x = element_blank())

ggsave(paste0(figure_dir, "growth_rates.png"), p, width = 9.52, height = 4.2, units = "cm", dpi = 300)

#median growth rates
medians <- sapply(growth_rates,median)
#checking how this compares with deterministic exp skyline growth
medians
#handwavy way of getting to population size under deterministic exp skyline growth
exp(medians[1]*8.333)*exp(medians[2]*8.333)*exp(medians[3]*8.333)
#another way: take the mean between all tree periods:
mean_over_skyline <- mean(medians)
mean_over_skyline
#calculate population under simple exp growth under that:
exp(mean_over_skyline*25)

#########################################
##plotting growth rates as skyline plot##
#########################################

typewriter_mcmc <- as.mcmc(typewriter)
#plotting the growth rate
growth <- typewriter_mcmc[,paste0(rep("birthRate.",13),1:13)] - typewriter_mcmc[,paste0(rep("deathRate.",13),1:13)]

HPD <- HPDinterval(growth)

median <- as.numeric(sapply( data.frame(growth),median))

#deterministic population size
total <- 1
for(i in 1:12) {

  total <- total*exp(median[i]*2)

}
total <- total*exp(median[13])

up_bd <- as.numeric(HPD[,"upper"])
low_bd <- as.numeric(HPD[,"lower"])

median <- c(median,median[length(median)])
up_bd <- c(up_bd,up_bd[length(up_bd)])
low_bd <- c(low_bd,low_bd[length(low_bd)])

#formatting the timeline into dates:
timeline_format <- c(seq(0,25,by=2),25)

#creating a dataframe and formatting for step plot
data_growth <- data.frame(timeline_format,median,low_bd,up_bd)
colnames(data_growth) <- c("Date","Median","95% HPI lower","95% HPI upper")
df_growth <- pivot_longer(data_growth,c("Median","95% HPI lower","95% HPI upper"),names_to = "stat",values_to = "value")
type <- df_growth$stat
type[which(type == "95% HPI lower")] <- "95% HPI"
type[which(type == "95% HPI upper")] <- "95% HPI"
df_growth <- cbind(df_growth,type)

p_growth <- ggplot(df_growth, aes(x=Date, y = value, key = stat, linetype = type)) +
  geom_step(size=1) +
  scale_color_manual(values=c("#5CA17D","#5CA17D","#5CA17D")) +
  scale_linetype_manual(values=c("dotted","solid")) +
  xlab("Time") +
  ylab(parse(text = paste0('"Growth rate "', '(~d^-1)')))  +
  theme(text=element_text(size = text_size),
          axis.line = element_line(color="black", size = 0.5),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = c(0.7,0.7),
          legend.title = element_blank())
p_growth
ggsave(paste0(figure_dir, "growth_skyline.png"), p_growth , width = 9.52, height = 4.2, units = "cm", dpi = 300)


