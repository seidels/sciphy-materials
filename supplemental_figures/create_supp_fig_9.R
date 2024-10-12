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

## figure settings
figure_path = "plots/supp_fig_9.pdf"

## load up the packages we will need:  
library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)

#load the log file
## NOTE the *2*-combined.log indicates to use the log of analysis 2 from figure 3.
typewriter_file <- "../figure_3/inference_output/2-combined.log"
typewriter <- read.table(typewriter_file, header = T) 

#########################################
##plotting growth rates as skyline plot##
#########################################

typewriter_mcmc <- as.mcmc(typewriter)
#plotting the growth rate
growth <- typewriter_mcmc[,paste0(rep("birthRate.",13),1:13)] - 
  typewriter_mcmc[,paste0(rep("deathRate.",13),1:13)]

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

p_growth <- ggplot(df_growth, aes(x=Date, y = value, key = stat, linetype = type )) +
  geom_step(size=1) + 
  scale_color_manual(values=c("#5CA17D","#5CA17D","#5CA17D")) +
  scale_linetype_manual(values=c("dotted","solid")) + 
  xlab("Time") + 
  ylab(parse(text = paste0('"Growth rate 95% HPI "', '(~day^-1)')))  + 
  theme(text=element_text(size = 22),axis.line = element_line(color="black", size = 0.5),
          panel.border = element_blank(),
          panel.background = element_blank(),panel.grid.major.x = element_blank(),legend.position = c(0.7,0.7),legend.title = element_blank()) 

ggsave(figure_path, p_growth , width = 50, height = 15, units = "cm", dpi = 800)


