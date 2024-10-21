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
library(RColorBrewer)
library(viridis)

#load the combined log file

typewriter_file <- "../figure_5/inference_output/4-mGASv2-skyline-ou.10burnin.combined.log"
typewriter <- read.table(typewriter_file, header = T)

figure_dir = "plots/"

# -----------------------------
#  plot insertion probabilities
# -----------------------------
substrLeft <- function(x, n){
  substr(x, 1, n)
}

#extract insertion probabilities from the full dataframe
insert_probs <- typewriter[,startsWith(x = colnames(typewriter), prefix = "editProbabilities.")]
insert_file = "../figure_5/processed_data/integer_conversion_table.csv"
insert_map = read.csv(insert_file)
trinucleotides_names <- insert_map$edits[which(insert_map$edits != "None")]
trinucleotides_names = unname(sapply(trinucleotides_names, function(x){substrLeft(x, n = 3)}))

#reorder by median
names(insert_probs) <- trinucleotides_names
insert_probs <- insert_probs[order(-sapply(insert_probs, median))]

# create custom colour palette
first_9_colors <- brewer.pal(9, "Set1")

# Generate 33 similar shades using viridis
remaining_33_colors <- viridis(33, option = "C")

# Combine the two sets of colors
custom_palette <- c(first_9_colors, remaining_33_colors, "white")

insert_color_table = data.frame(inserts = c(colnames(insert_probs), "None"), colors = custom_palette)
write.csv(insert_color_table, file = "plots/insert_colour_table.csv", quote = F, row.names = F)

# add insert prob posterior summaries
insert_probs <- bind_cols(insert_probs, Prior=dirmult::rdirichlet(nrow(insert_probs), rep(1.5,19))[,2])
insert_probs_medians <- sapply(insert_probs, median)
insert_probs_low <- as.numeric(sapply(insert_probs,function(x) {quantile(x, 0.025)}))
insert_probs_up <-  as.numeric(sapply(insert_probs,function(x) {quantile(x, 0.975)}))
datafra <- data_frame(name=c(trinucleotides_names,"Prior"), median=insert_probs_medians,low=insert_probs_low,up=insert_probs_up)


datafra$name <- factor(datafra$name, levels = datafra$name)

# Print the custom palette to check
#print(custom_palette)

p_inserts <-  ggplot(datafra) +
  geom_bar(aes(y=name,x=median, fill=name),stat="identity",colour="black") +

  theme_classic() +
  xlab("Posterior insert probability") + ylab("Insert") +
  theme(legend.position = "none" ) +
  scale_fill_manual(values= custom_palette) +
  #theme(plot.title = element_text(hjust=0.5)) +
  geom_errorbar(aes(y=name, xmin=low, xmax=up), width=.2,
                position=position_dodge(.9)) +
  #ggtitle("Insert probabilities") +

  coord_cartesian(xlim=c(0,0.35),expand = FALSE) +
  theme(
        axis.text = element_text(size = (text_size-2)),
        axis.title.y = element_text(size = text_size))

p_inserts


ggsave(paste0(figure_dir, "insert_probs.png"), p_inserts, width = 6, height = 11, units = "cm", dpi = 300)
