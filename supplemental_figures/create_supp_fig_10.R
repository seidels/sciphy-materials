## Script name: create_supp_fig_16
##
## Purpose of script: Plot a figure that summarises how SciPhy compares with UPGMA in terms of growth rate estimation
##


## set output directory for the plot
figure_path = "plots/supp_fig_16.pdf"
# set figure settings
text_size=10
label_size = 12

## load up the packages we will need:  (uncomment as required)
library(tidyverse)
library(data.table)
library(ape)
library(TreeDist)
library(gdata)
library(reshape2)
library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)
library(pammtools)
library(reshape2)


# ---------------------------------------------------------
#  plot the constant growth rates inferred on fixed trees and co-estimated with SciPhy trees
# ---------------------------------------------------------

## Add the UPGMA based estimates
typewriter_file <- "../figure_4/inference_output/fixed_tree/combined_fixed_UPGMA_rho_sampling.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
HDInterval::hdi(bd_rates)
bd_rates$tree = "UPGMA ordered"

bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = bd_rates

## Add UPGMA scaled based estimates
typewriter_file <- "../figure_4/inference_output/fixed_tree/combined_1000_UPGMA_infer_rho_sampling_median_height.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA ordered + \n root scaling"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the main SciPhy analysis with co-estimated rho sampling
typewriter_file <- "../figure_3/inference_output/combined_clockPerTarget_sampling_DataSet1.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "SciPhy"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the SciPhy analysis with fixed rho sampling
typewriter_file <- "../figure_3/inference_output/combined.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "SciPhy + \n Fixed sampling"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the analysis based on the UPGMA tree rescaled with SciPhy
typewriter_file <- "../figure_4/inference_output/fixed_tree/combined_1000_UPGMA_medianPosteriorHeight_estimateBranchLengths_infer_rho_sampling.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA ordered + \n SciPhy scaling"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)


## Add the Prior
prior = data.frame(growthRate=rlnorm(n = 1000, meanlog = 0.1, sdlog = 1) -  rlnorm(1000, meanlog = -0.4, sdlog = 1))
prior$tree = "Prior"
prior =melt(prior, id.vars = "tree")

growth_rates = rbind(growth_rates, prior)

# Define order
growth_rates$tree = factor(growth_rates$tree, levels=c("Prior", "SciPhy","SciPhy + \n Fixed sampling", "UPGMA ordered","UPGMA ordered + \n SciPhy scaling", "UPGMA ordered + \n root scaling"))
# Define colors
cols <- c("Prior" = "white","SciPhy + \n Fixed sampling"="#B7CC62","SciPhy" = "#56996E", "UPGMA ordered" = "red","UPGMA ordered + \n root scaling" = "#E07E5E","SciPhy MCC" = "grey","UPGMA ordered + \n SciPhy scaling" = "yellow","UPGMA ordered + \n root scaling" = "pink","SciPhy MCC"="grey")

p_growth <- ggplot(growth_rates, aes(x=tree, y=value, fill = tree)) +
  theme_bw() +
  geom_violin(draw_quantiles = c(0.5))+  
  scale_fill_manual(values = cols) + 
  scale_y_continuous(breaks=c(0.2,0.6,1.0,1.4),limits=c(0.2, 1.6)) + 
  xlab("Tree used for inference")+
  ylab(expression("Growth rate [" * d^-1 * "]")) + theme(legend.title=element_blank(),axis.text.x = element_text(size=7,color = "black") )

ggsave(figure_path,p_growth + theme(legend.position = "none") , width = 15, height = 7.14, units = "cm", dpi = 300)

