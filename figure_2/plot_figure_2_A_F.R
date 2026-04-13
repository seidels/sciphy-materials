## ---------------------------
##
## Script name: get plot for all inferred parameters in the full validation
##
## Purpose of script: plot well-calibrated simulations
##
## Author: Antoine Zwaans, Sophie Seidel
##
## Date Created: 2024-07-05
##
## Copyright (c) Antoine Zwaans, Sophie Seidel, 2024
## Email: azwaans@ethz.ch, sophie.seidel@posteo.de
##
## ---------------------------

#Guidelines Setup ---
text_size_max <- 7
text_size_min <- 5
plot_text_size <- 7 
axis_line_weight <- 0.2

#  Load Packages ---
require(data.table)
library(tidyverse)
library(ape)
library(treebalance)
library(TreeSim)
library(rstatix)
library(ggpubr)
library(phylobase)
library(tracerer)
library(HDInterval)
library(ggplot2)
library(stringr)
library(phytools)
library(cowplot)
library(TreeDist)

#  Read-in validation results ---
B1_index_inference <- read.csv(file = "inference_logs/summary/B1_index_inference.csv")
insert_rate_inference <- read.csv(file = "inference_logs/summary/insert_rate_inference.csv")
clock_rate_inference <- read.csv(file = "inference_logs/summary/clock_rate_inference.csv")
tree_height_inference <- read.csv(file = "inference_logs/summary/tree_height_inference.csv")
tree_length_inference <- read.csv(file = "inference_logs/summary/tree_length_inference.csv")

nr_converged_chains <- 100

#  Coverage and correlation statistics ---
coverages_per_insert <- c()
for(i in 1:13) {
  coverage <- sum(insert_rate_inference[which(insert_rate_inference$insertRate == i),"recovered"])/nr_converged_chains
  coverages_per_insert <- c(coverages_per_insert, coverage)
}

coverage_clock <- sum(clock_rate_inference[,"recovered"])/100
coverage_height <- sum(tree_height_inference[,"recovered"])/100
coverage_length <- sum(tree_length_inference[,"recovered"])/100
coverage_B1 <- sum(B1_index_inference[,"recovered"])/100


correlation_clock <- cor.test(clock_rate_inference$median, clock_rate_inference$true_value, method = "pearson")
correlation_tree_height <- cor.test(tree_height_inference$median, tree_height_inference$true_value, method = "pearson")
correlation_tree_length <- cor.test(tree_length_inference$median, tree_length_inference$true_value, method = "pearson")
correlation_B1 <- cor.test(B1_index_inference$median, B1_index_inference$true_value, method = "pearson")

# Theme set
theme_set(
  theme_classic(base_size = plot_text_size, base_family = "") +
    theme(
      plot.title = element_text(hjust = 0.5, size = text_size_max, face = "bold"),
      axis.title = element_text(size = text_size_max),
      axis.text = element_text(size = text_size_min),
      axis.line = element_line(linewidth = axis_line_weight),
      axis.ticks = element_line(linewidth = axis_line_weight),
      legend.position = "none"
    )
)

cols_recovered <- c("TRUE" = "black", "FALSE" = "darkgrey")

# Panels A to E: Generate parameter Validation Plots ---
make_val_plot <- function(df, title, x_lim, y_lim, y_lab = "Estimated median\n& posterior interval") {
  ggplot(df, aes(x=true_value, y=median, color=as.logical(recovered))) +
    geom_point(size = 0.2) + 
    geom_errorbar(aes(ymin = hpd_lower, ymax = hpd_upper), alpha = 0.4, linewidth = 0.2) +
    geom_abline(slope = 1, col = "darkgreen", linewidth = 0.3) +
    scale_color_manual(values = cols_recovered) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(title = title, x = "True value", y = y_lab)
}

insert_probs_plot <- make_val_plot(insert_rate_inference, "Insertion probabilities", c(0, 0.5), c(0, 0.5))
clock_rate_plot   <- make_val_plot(clock_rate_inference, "Editing rate", c(0, 0.5), c(0, 0.5))
tree_height_plot  <- make_val_plot(tree_height_inference, "Tree height", c(15, 30), c(15, 30))
tree_length_plot  <- make_val_plot(tree_length_inference, "Tree length", c(0, 10000), c(0, 10000))
B1_index_plot     <- make_val_plot(B1_index_inference, "Tree balance", c(0, 220), c(0, 220))

# Panel F - Tree Topology Analysis 
distance_SCIPHY_truth_PI <- c()
distance_truth_random_PI <- c()

for(SEED in 1:100) {
  clean_tree <- function(path) {
    tr <- readLines(path)
    tr <- str_remove_all(tr, "\\[.*?\\]") 
    tr <- paste0(tr, ";")
    ape::read.tree(text = tr)
  }
  
  CCD_tree <- ape::read.nexus(paste0("inference_logs/CCD/CCD_tree.", SEED, ".txt"))
  true_tree <- clean_tree(paste0("simulated_data/simulate_alignment_and_tree.", SEED, ".newick"))
  
  random_tree <- TreeSim::sim.bd.taxa.age(length(true_tree$tip.label), 1, 0.8, 0.2, 0.00003, 25)
  random_tree[[1]]$tip.label <- str_remove(random_tree[[1]]$tip.label, "t")
  
  distance_SCIPHY_truth_PI <- c(distance_SCIPHY_truth_PI, TreeDist::PhylogeneticInfoDistance(true_tree, CCD_tree, normalize = TRUE))
  distance_truth_random_PI <- c(distance_truth_random_PI, TreeDist::PhylogeneticInfoDistance(random_tree[[1]], true_tree, normalize = TRUE))
}

pi_distances_to_truth <- data.frame(seed=1:100, SciPhy_CCD=distance_SCIPHY_truth_PI, Random_BD_tree=distance_truth_random_PI) %>%
  pivot_longer(!seed, names_to = "reference", values_to = "distance")

# Clean the data by removing incomplete pairs
pi_distances_filtered <- pi_distances_to_truth %>%
  group_by(seed) %>%
  # Keep only seeds where NO distance value is NaN
  filter(!any(is.na(distance))) %>%
  ungroup()


stat.test <- pi_distances_filtered %>%
  pairwise_t_test(distance ~ reference, paired = TRUE, p.adjust.method = "bonferroni", detailed = TRUE) %>%
  add_xy_position(x = "reference")


effect_size <- pi_distances_filtered %>%
  cohens_d(distance ~ reference, paired = TRUE)


print(stat.test)
print(effect_size)

#  Boxplot for Topology ---
cols_topo <- c("SciPhy_CCD" = "#2e7d32cc", "Random_BD_tree" = "#005bf2cc")

bxp_CCD_PI <- ggboxplot(
  pi_distances_to_truth, x = "reference", y = "distance",
  fill = "reference", linewidth = 0.2, outlier.size = 0.3
) + 
  scale_fill_manual(values = cols_topo, labels = c("SciPhy CCD", "Random BD tree")) +
  scale_x_discrete(labels = c("SciPhy CCD", "Random BD tree")) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  coord_cartesian(ylim = c(0.0, 1.3)) + 
  labs(title = "Tree topology", x = "Tree reconstruction method", y = "Normalized PI distance \n to true tree")

# Updated Boxplot Chunk
bxp_CCD_PI <- ggboxplot(
  pi_distances_to_truth, 
  x = "reference", 
  y = "distance",
  fill = "reference", 
  linewidth = 0.2, 
  outlier.size = 0.2 
) + 
  scale_fill_manual(
    values = cols_topo, 
    # Ensure these are in the order as the 'reference' column
    labels = c("Random_BD_tree" = "Random BD tree", "SciPhy_CCD" = "SciPhy CCD")
  ) +
  scale_x_discrete(
    labels = c("Random_BD_tree" = "Random BD tree", "SciPhy_CCD" = "SciPhy CCD")
  ) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  coord_cartesian(ylim = c(0.0, 1.3)) + 
  labs(
    title = "Tree topology", 
    x = "Tree reconstruction method", 
    y = "Normalized PI distance \n to true tree"
  )

# match axis thickness of other panels
paired_test_CCD_PI <- bxp_CCD_PI + 
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif", 
    step.increase = 0.1, 
    y.position = 1.1,
    size = 2 
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = text_size_max, face = "bold"),
    axis.title = element_text(size = text_size_max),
    axis.text = element_text(size = text_size_min),
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(linewidth = axis_line_weight),
    legend.position = "none"
  )
# Assembly in a grid
full_figure <- cowplot::plot_grid(
  insert_probs_plot, clock_rate_plot, 
  tree_height_plot, tree_length_plot, 
  B1_index_plot, paired_test_CCD_PI, 
  ncol = 2, 
  labels = "AUTO",
  label_size = text_size_max
)

ggsave(
  filename = "plots/figure_2.pdf",
  plot = full_figure,
  width = 180, 
  height = 185, 
  units = "mm", 
  device = "pdf", 
  dpi = 300
)
