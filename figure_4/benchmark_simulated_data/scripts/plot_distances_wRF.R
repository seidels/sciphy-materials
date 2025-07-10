## ---------------------------
##
## Purpose of script: Plot benchmark
##
## Author: Sophie Seidel
##
## Date Created: 2025-06-15
##
## Copyright (c) Sophie Seidel, 2025
## Email: soseidel@uw.edu
##

# libs
library(ggplot2)
library(tidyverse)
library(wesanderson)
palette = c("#2E7D32", wes_palette("Zissou1", n = "5")[c(3, 1, 4, 5)])

# read data
dat_file = "benchmark/distances_across_methods.RDS"

dat = readRDS(dat_file)
columns = colnames(dat)

# select columns to plot
columns_to_transform = columns[startsWith(x = columns, prefix = "wRF")]
columns_to_transform = setdiff(columns_to_transform, c("wRF_distance_sciphy_mcc", "wRF_distance_tidetree_mcc"))

# define plotting order
method_order <- c(
  "wRF_distance_sciphy_ccd",
  "wRF_distance_upgma_ordered",
  "wRF_distance_tidetree_ccd",
  "wRF_distance_upgma",
  "wRF_distance_random"
)

dat_long = dat %>% 
  pivot_longer(cols = columns_to_transform) %>%
  mutate(name = factor(name, levels = method_order))



p = ggplot(dat_long, aes(x = name, y = value, fill = name)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = palette)+
  scale_x_discrete(labels = c(
    "wRF_distance_sciphy_ccd" = "SciPhy CCD",
    "wRF_distance_tidetree_ccd" = "TiDeTree CCD",
    "wRF_distance_upgma" = "UPGMA unordered",
    "wRF_distance_upgma_ordered" = "UPGMA ordered",
    "wRF_distance_random" = "Random BD tree"
  ))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "",
    y = "Normalised wRF distance",
    fill = "Method"
  ) +
  theme_classic() +
  theme(text = element_text(size = 12), 
        legend.position = "none")

p
ggsave(filename = "benchmark/plots/method_comparison_wRF.pdf", plot = p, dpi = 300, 
       width = 18, height = 8, units = "cm")

svg(filename = "benchmark/plots/method_comparison_wRF.svg", 
    width = 7.1, height = 3.1)
p
dev.off()
ggsave(filename = "benchmark/plots/method_comparison_wRF.jpg", plot = p, dpi = 300, 
       width = 18, height = 10, units = "cm")
