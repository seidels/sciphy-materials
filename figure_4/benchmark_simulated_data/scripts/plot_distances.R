
library(ggplot2)
library(tidyverse)
library(wesanderson)
palette = c("#2E7D32", wes_palette("Zissou1", n = "5")[c(1,3:4, 5)])

# read data
dat_file = "benchmark/distances_across_methods.RDS"

dat = readRDS(dat_file)
columns = colnames(dat)

# select columns to plot
columns_to_transform = columns[startsWith(x = columns, prefix = "RF")]
columns_to_transform = setdiff(columns_to_transform, c("RF_distance_sciphy_mcc", "RF_distance_tidetree_mcc"))

# define plotting order
method_order <- c(
  "RF_distance_sciphy_ccd",
  "RF_distance_tidetree_ccd",
  "RF_distance_upgma_ordered",
  "RF_distance_upgma",
  "RF_distance_random"
)


dat_long = dat %>% 
  pivot_longer(cols = columns_to_transform) %>%
  mutate(name = factor(name, levels = method_order))



p = ggplot(dat_long, aes(x = name, y = value, fill = name)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = palette, labels = c(
                                                 "RF_distance_sciphy_mcc" = "SciPhy MCC",
                                                 "RF_distance_tidetree" = "TiDeTree MCC",
                                                 "RF_distance_upgma" = "UPGMA unordered",
                                                 "RF_distance_upgma_ordered" = "UPGMA ordered",
                                                 "RF_distance_random" = "Random tree"))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Tree Size Category",
    y = "RF distance",
    fill = "Method"
  ) +
  theme_bw() +
  theme(text = element_text(size = 12), 
        legend.position = "top")

p
ggsave(filename = "benchmark/plots/method_comparison.pdf", plot = p, dpi = 300, 
       width = 20, height = 10, units = "cm")
ggsave(filename = "benchmark/plots/method_comparison.jpg", plot = p, dpi = 300, 
       width = 20, height = 10, units = "cm")
