## ---------------------------
##
## Script name: create_supp_fig_15
##
## Purpose of script: Plot estimates from SciPhy on HEK293 cell culture data, skyline growth + clock per target, with color annotations matching the tree.
##
## Author: Antoine Zwaans
##
## Date Created: 2023-10-13
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch

## figure settings
figure_path = "plots/supp_fig_15.pdf"

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

text_size_max <- 10
text_size_min <- 8
line_thickness <- 0.5
growth_color <- "#B7CC62" 

p_growth <- ggplot() +  
  geom_stepribbon(data = data_growth, 
                  aes(x = Date, ymin = `95% HPI lower`, ymax = `95% HPI upper`, fill = "SciPhy"), 
                  alpha = 0.2) + 
  
  geom_step(data = subset(df_growth, stat == "Median"), 
            aes(x = Date, y = value, color = "SciPhy"), 
            linewidth = 1) + 
  
  scale_fill_manual(values = c("SciPhy" = growth_color)) +
  scale_color_manual(values = c("SciPhy" = growth_color)) +
  
  labs(y = expression("Growth rate [" * d^-1 * "]"), 
       x = "Time [d]") +
  
  theme_classic(base_size = text_size_min) + 
  theme(
    legend.position = "none", 
    legend.title = element_blank(), 
    legend.background = element_blank(),
    axis.title = element_text(size = text_size_max),
    axis.text = element_text(size = 7),
    axis.line = element_line(linewidth = line_thickness), 
    axis.ticks = element_line(linewidth = line_thickness),
    panel.grid.major.x = element_blank()
  )


ggsave(figure_path,p_growth + theme(legend.position = "none") , width = 15 , height = 7.14, units = "cm", dpi = 300)
