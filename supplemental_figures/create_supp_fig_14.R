## ---------------------------
##
## Script name: create_supp_fig_20
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


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(dirmult)
library(grid)

text_size_max <- 10
text_size_min <- 8
axis_line_weight <- 0.4
figure_dir <- "plots/"

typewriter_file <- "../figure_5/inference_output/4-mGASv2-skyline-ou.log"
typewriter <- read.table(typewriter_file, header = TRUE)

substrLeft <- function(x, n) substr(x, 1, n)

insert_probs <- typewriter[, startsWith(x = colnames(typewriter), prefix = "editProbabilities.")]
insert_file <- "../figure_5/processed_data/integer_conversion_table.csv"
insert_map <- read.csv(insert_file)

trinucleotides_names <- insert_map$edits[which(insert_map$edits != "None")]
trinucleotides_names <- unname(sapply(trinucleotides_names, function(x) substrLeft(x, n = 3)))
names(insert_probs) <- trinucleotides_names

# add prior
set.seed(42) 
n_samples <- nrow(insert_probs)
prior_samples <- rdirichlet(n_samples, rep(1.5, 19))[, 2]
insert_probs$Prior <- prior_samples

# sort by median by median 
insert_medians <- sapply(insert_probs %>% select(-Prior), median)
sorted_insert_names <- names(sort(insert_medians, decreasing = TRUE))

# combine with Prior at the very end
ordered_names <- c(sorted_insert_names, "Prior")

# setup Colors
first_9_colors <- brewer.pal(9, "Set1")
remaining_33_colors <- viridis(33, option = "C")
custom_palette <- c(first_9_colors, remaining_33_colors, "white")

# map colors to ordered names, ensuring Prior is assigned its specific color
plot_colors <- setNames(custom_palette[1:length(ordered_names)], ordered_names)
plot_colors["Prior"] <- "#E1E1F7"


datafra_long <- insert_probs %>%
  pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
  mutate(name = factor(name, levels = ordered_names))

# generate Violin Plot ---
p_inserts <- ggplot(datafra_long, aes(x = name, y = value, fill = name)) +
  geom_violin(scale = "width", draw_quantiles = c(0.5), linewidth = 0.2, adjust = 2.2) +
  scale_fill_manual(values = plot_colors) +
  
  coord_cartesian(ylim = c(0, 0.35)) +
  scale_y_continuous(
    breaks = seq(0, 0.35, 0.05), 
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(y = "Insertion probability",x="Estimates per trinucleotide insert") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.1),
    axis.ticks = element_line(linewidth = axis_line_weight),
    axis.title.y = element_text(size = text_size_max),
    axis.text.y = element_text(size = text_size_min),
    axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size = text_size_min),
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  ) 


p_inserts

ggsave(
  filename = "plots/supp_fig_20.pdf", 
  plot = p_inserts, 
  width = 15, 
  height = 14, 
  units = "cm", 
  device = cairo_pdf 
)
