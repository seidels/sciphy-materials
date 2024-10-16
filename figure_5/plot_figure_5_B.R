## ---------------------------
##
## Script name: plot_inference_results_annotated
##
## Purpose of script: Plot estimates from SciPhy on HEK293 cell culture data, clock per target, with color annotations matching the tree.
##
## Author: Antoine Zwaans & Sophie Seidel
##
## Date Created: 2023-10-13
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##


# set figure settings
text_size=10
label_size = 12
## load up the packages we will need:

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)
library(magick)

#load the combined log file

typewriter_file <- "inference_output/4-mGASv2-skyline-ou.10burnin.combined.log"
typewriter <- read.table(typewriter_file, header = T)

figure_dir = "plots/"

#extract the clock rate from the tract
clock_rate <- typewriter[, startsWith(colnames(typewriter), "clockRate_")]

#rename by the targetBC
order_of_tbcs_in_alignment =c("AGGCTAATTCCC", "TAACGAAGATTT", "AACTGAATGTTT",
                              "GTATAAAGTTTG", "TTGATAACGTGA", "GTTGAAAGGTGA", 
                              "GTACAAAATAGT", "TTCACAAATTTA") 

tape_to_int_map = read.csv(file = "tape_seq_identifier_map.csv")
numeric_names = unname(sapply(order_of_tbcs_in_alignment, function(x){
  which(tape_to_int_map$tape_identifier == x)
}))

names(clock_rate) <- numeric_names

#add a prior column
clock_rate <- bind_cols(clock_rate,Prior=rlnorm(nrow(clock_rate), meanlog = -2, sdlog = 0.5))
clock_rate_long <- pivot_longer(clock_rate,seq(1,ncol(clock_rate)))

#order columns
#clock_rate_long <- mutate(clock_rate_long, name = fct_relevel(name,ordered_names,"Prior"))

p_clock_pos <- ggplot(clock_rate_long, aes(x=name,value,fill=name)) +
  theme_classic() +
  geom_violin(draw_quantiles =  c(0.5), width=2) +
  ylab(expression("Posterior edit rate [" * d^-1 * "]")) +
  xlab("Tapes") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0, 0.8))+
  theme(legend.position = "none",
    axis.title = element_text(size = text_size),
        axis.text.y = element_text(margin = margin(r = 0, l=-0.5)),
        axis.text = element_text(size = (text_size-2))
    ) 

p_clock_pos

#ggsave(paste0(pic_dir,"clock_rate.png"), p_clock_pos, width = 4.76, height = 7.14, units = "cm", dpi = 300)
ggsave(paste0(figure_dir,"edit_rate.png"), p_clock_pos, width = 7, height = 4.7, units = "cm", dpi = 300)
  