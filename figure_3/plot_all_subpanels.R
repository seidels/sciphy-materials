## ---------------------------
##
## Script name: plot_inference_results_annotated
##
## Purpose of script: Plot estimates from SciPhy on HEK293 cell culture data, 
## clock per target, with descending median ordering.
##
## Author: Antoine Zwaans & Sophie Seidel
##
## Date Created: 2023-10-13
## Updated for Compliance: 2026-04-03
##

# --- 1. Global Constraints ---
text_size_max <- 7
text_size_min <- 5
axis_line_weight <- 0.2
figure_dir = "plots/"

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(scales)
library(grid)
library(ggtree)
library(treeio)
library(cowplot)

# Load data
typewriter_file <- "inference_output/3-combined.log"
typewriter <- read.table(typewriter_file, header = T)

# -----------------------------
#  Plot Insertion Probabilities
# -----------------------------

insert_probs <- typewriter[,c(7:25)]
trinucleotides_names <- c("CAT","CCG","GCC","ATA","GAT","ACG","ACA","TCG","TAT","GCT","CTA","TGT","AGA","TAA","CAG","TAG","GAG","ACC","GCG")
names(insert_probs) <- trinucleotides_names

# Order by descending median
insert_probs <- insert_probs[order(-sapply(insert_probs, median))]
ordered_trinucs <- names(insert_probs)

insert_probs <- bind_cols(insert_probs, Prior = rdirichlet(nrow(insert_probs), rep(1.5, 19))[,2])

# Pivot to long format for geom_violin
datafra_long <- insert_probs %>%
  pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
  mutate(name = factor(name, levels = names(insert_probs)))

n_categories <- nlevels(datafra_long$name)
estimate_span <- n_categories - 1
colors <- c(hue_pal()(19), "#E1E1F7")
names(colors) <- c(sort(trinucleotides_names), "Prior")

p_inserts <- ggplot(datafra_long, aes(x = name, y = value, fill = name)) +
  geom_violin(scale = "width", draw_quantiles = c(0.5), linewidth = 0.2, trim = TRUE, bounds = c(0, Inf)) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 0.2), clip = "off") +
  annotation_custom(grid::textGrob(label = "Estimates per insert", x = grid::unit((1 + estimate_span)/2, "native"),
                                   y = grid::unit(-55, "points"), gp = grid::gpar(fontsize = text_size_max))) +
  labs(y = "Insertion probability") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(linewidth = axis_line_weight),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = text_size_max),
    axis.text.y = element_text(size = text_size_min),
    axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1.0, size = text_size_min),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5)
  )

# --------------------
#  Plot Clock Rates
# --------------------

clock_rate_data <- typewriter[,c(29:41)]
names(clock_rate_data) <- c("ATGGTAAG","ATTTATAT","ATTTGGTT", "GCAGGGTG","GTAAAGAT", 
                            "TAGATTTT","TGCGATTT", "TGGACGAC","TGGTTTTG", "TTAGATTG",
                            "TTGAGGTG","TTTCGTGA","TTCACGTA")

ordered_clocks <- names(clock_rate_data)[order(-sapply(clock_rate_data, median))]
clock_rate_data <- bind_cols(clock_rate_data, Prior = rlnorm(nrow(clock_rate_data), -2, 0.5))
clock_rate_long <- pivot_longer(clock_rate_data, everything()) %>%
  mutate(name = fct_relevel(name, ordered_clocks, "Prior"))

segments <- data.frame(xstart = c(1, 2, 5), xstop = c(1, 4, 13), ystart = c(0.35, 0.25, 0.21), ystop = c(0.35, 0.25, 0.21))

p_clock_pos <- ggplot(clock_rate_long, aes(x = name, value, fill = name)) +
  geom_violin(draw_quantiles = c(0.5), linewidth = 0.2, trim = TRUE) +
  scale_fill_manual(values = c(c("#0A0A82", rep("#404085", 3), rep("#7C7CA3", 9)), "#E1E1F7")) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0.1, 0.37), clip = "on") + 
  labs(y = expression("Editing rate [" * d^-1 * "]")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(linewidth = axis_line_weight),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = text_size_max),
    axis.text.y = element_text(size = text_size_min),
    axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1, size = text_size_min),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5)
  ) +
  annotate("text", x = c(1, 3, 9), y = c(0.36, 0.27, 0.24), label = c("2X", "4X", "5X"), size = text_size_min * 0.3527) +
  geom_segment(data = segments, aes(x = xstart, xend = xstop, y = ystart, yend = ystop), inherit.aes = FALSE, linewidth = 0.2)

# --------------------
#  Growth Rate & Sampling Proportion
# --------------------

bd_rates <- typewriter[,"birthRate"] - typewriter[,"deathRate"]
bd_rates <- bind_cols(Estimate = bd_rates, Prior = rlnorm(length(bd_rates), -0.6, 1) - rlnorm(length(bd_rates), -2, 1))
bd_rates_long <- pivot_longer(bd_rates, everything()) %>% mutate(name = fct_relevel(name, "Estimate", "Prior"))

p_growth <- ggplot(bd_rates_long, aes(x = name, value, fill = name)) +
  geom_violin(draw_quantiles = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = c("white", "#E1E1F7")) +
  coord_cartesian(ylim = c(-0.1, 0.9)) +
  labs(y = expression("Growth rate [" * d^-1 * "]")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(linewidth = axis_line_weight),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = text_size_max),
    axis.text = element_text(size = text_size_min)
  )

samp_prop <- typewriter[,"samplingProportion"] 
samp_prop <- bind_cols(Estimate = samp_prop, Prior = rlnorm(length(samp_prop), -7.4, 1.2))
samp_prop_long <- pivot_longer(samp_prop, everything()) %>% mutate(name = fct_relevel(name, "Estimate", "Prior"))

p_samp <- ggplot(samp_prop_long, aes(x = name, y = value, fill = name)) +
  geom_violin(draw_quantiles = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = c("white", "#E1E1F7")) +
  coord_cartesian(ylim = c(0.0, 0.0005)) +
  labs(y = "Sampling proportion") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(linewidth = axis_line_weight),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = text_size_max),
    axis.text = element_text(size = text_size_min)
  )

# --------------------
#  Phylogenetic Trees
# --------------------

tree = treeio::read.beast(file = tree_file)
stem_length = 25 - max(nodeHeights(tree@phylo))

p <- ggtree(tree, root.position = stem_length, size=0.05) +
  geom_hilight(node=1309, fill="steelblue", alpha=0.5) +
  geom_rootedge(rootedge = stem_length, size=0.05) + 
  theme_tree2() +
  theme(
    axis.line.x = element_line(linewidth = axis_line_weight),
    axis.ticks.x = element_line(linewidth = axis_line_weight),
    axis.text.x = element_text(size = text_size_min),
    axis.title.x = element_text(size = text_size_max)
  ) +
  xlab("Time [d]")

# 2. Extract plot data to find the visual range
d <- p$data
# Get all offspring of node 1309
clade_nodes <- tidytree::offspring(tree, 1309)

# 3. Find the tips with the extreme Y-coordinates within that clade
clade_tips <- d[d$node %in% clade_nodes & d$isTip, ]
taxa1_node <- clade_tips$node[which.min(clade_tips$y)]
taxa2_node <- clade_tips$node[which.max(clade_tips$y)]

p <- p + geom_strip(
  taxa1 = taxa1_node, 
  taxa2 = taxa2_node, 
  color = "steelblue", 
  barsize = 1,
  offset = 0.5,     # Adjust to move the bracket away from the tips
  extend = 0.2      # Adjust to make the bracket slightly taller/shorter
)

# Clade Visualization with Heatmap
sub = tree_subset(tree = tree, node = 1309, group_node = T, root_edge = T, levels_back = F)
subtree_stem_length = 25 - max(nodeHeights(get.tree(sub)))

sub_tree_plot = ggtree(sub, root.position = subtree_stem_length) + 
  theme_tree2() + 
  geom_rootedge(rootedge = subtree_stem_length) +  
  geom_range('height_0.95_HPD', color='grey', size=1, alpha=.4) +
  geom_nodepoint(aes(size=posterior)) + 
  scale_size(range = c(0.1, 1.5), name = "Posterior support") +
  theme(
    axis.line.x = element_line(linewidth = axis_line_weight),
    axis.ticks.x = element_line(linewidth = axis_line_weight),
    axis.text.x = element_text(size = text_size_min),
    axis.title.x = element_text(size = text_size_max)
  )

# Heatmap loop logic
for (i in 1:length(targetBCs)){
  labels <- c(targetBCs[i], rep("", ncol(targetbc_edits_list[[i]])-1))
  sub_tree_plot = gheatmap(sub_tree_plot, targetbc_edits_list[[i]], 
                           colnames_position = "bottom", width = 0.04, offset = 1.1 * (i-1),
                           colnames_offset_y = -10,custom_column_labels = labels, font.size = (text_size_min / 3), 
                           colnames_angle = 70) + coord_cartesian(clip="off") 
}

final_tree <- sub_tree_plot +
  theme(
    legend.position = "top", 
    legend.box = "vertical", # Forces separate lines
    legend.margin = margin(b = -10),
    legend.text = element_text(size = text_size_min),
    legend.title = element_text(size = text_size_min),
    legend.key.size = unit(0.2, 'cm'),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5)
  ) +
  guides(
    fill = guide_legend(nrow = 2, title = "Trinucleotide inserts"), 
    size = guide_legend(nrow = 1)
  ) +
  xlab("Time [d]")

# --------------------
#  Combined Figure Export
# --------------------

blank_p <- ggplot() + theme_void()


# 1. First nesting level (A and B)
blank_and_clock <- plot_grid(
  blank_p, p_clock_pos, 
  nrow = 2, 
  rel_heights = c(1, 1.5), 
  labels = c("A", "B"),
  label_size = 7  # Forces 7pt labels
)

# 2. Second nesting level (C)
clock_and_insert_probs <- plot_grid(
  blank_and_clock, p_inserts, 
  ncol = 2, 
  labels = c("", "C"),
  label_size = 7
)

# 3. Third nesting level (D, E, F)
growth_and_sampling <- plot_grid(
  p_growth, p_samp, p, 
  ncol = 3, 
  rel_widths = c(1, 1, 2), 
  labels = c("D", "E", "F"),
  label_size = 7
)

# 4. Final assembly (G)
all_plots <- plot_grid(
  clock_and_insert_probs, 
  growth_and_sampling, 
  final_tree, 
  nrow = 3, 
  rel_heights = c(3, 1.5, 4), 
  labels = c("", "", "G"),
  label_size = 7
)

# Save with requested dimensions (185 mm x 180 mm)
ggsave(paste0(figure_dir, "all_plots_final_GUIDELINES_ADJUSTED.pdf"), all_plots, width = 180, height = 185, units = "mm",)

