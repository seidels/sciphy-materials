## ---------------------------
##
## Script name: create_supp_fig_6
##
## Purpose of script: Plot inference results for the analysis where only 2 sites were included
##
## Author: Sophie Seidel
##
## Date Created: 2024-03-18
##
## Copyright (c) Sophie Seidel, 2024
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


## set working directory where the log files are

data_dir = "supp_fig_6_data/"

## figure settings
text_size = 11
figure_path = "./plots/supp_fig_6.pdf"

## load up the packages we will need:

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)

#load the combined log file

log_file_2_sites <- "supp_fig_2_data/analysis_2sites_DataSet1.combined.burnin20.log"
dat_2_sites <- read.table(log_file_2_sites, header = T)

log_file_all_sites = "../figure_3/inference_output/1-combined.log"
dat_all_sites = read.table(log_file_all_sites, header = T)
dat_all_sites = dat_seq_data

# --------------------
#  plot the clock rates
# --------------------

#extract the clock rate from the tract
clock_rate_all_sites = dat_all_sites[, startsWith(x = colnames(dat_all_sites), prefix = "clock")]
clock_rate_2_sites <- dat_2_sites[, startsWith(x = colnames(dat_2_sites), prefix = "clock")]

#rename by the targetBC
names(clock_rate_2_sites) = names(clock_rate_all_sites) <- c("ATGGTAAG","ATTTATAT",
                       "ATTTGGTT", "GCAGGGTG",
                       "GTAAAGAT", "TAGATTTT",
                       "TGCGATTT", "TGGACGAC",
                       "TGGTTTTG", "TTAGATTG",
                       "TTGAGGTG","TTTCGTGA",
                       "TTCACGTA")


#reorder by median
ord_all = order(-sapply(clock_rate_all_sites, median))
clock_rate_all_sites <- clock_rate_all_sites[ord_all]

ord_2 = order(-sapply(clock_rate_2_sites, median))
clock_rate_2_sites <- clock_rate_2_sites[ord_2]

ordered_names_all <- names(clock_rate_all_sites)
ordered_names_2 <- names(clock_rate_2_sites)

#add a prior column
#clock_rate <- cbind(clock_rate,Prior=rlnorm(nrow(clock_rate), meanlog = -2, sdlog = 0.5))
clock_rate_all_sites_long <- pivot_longer(clock_rate_all_sites, seq(1,ncol(clock_rate_all_sites)))
clock_rate_all_sites_long$data = "all sites"
clock_rate_2_sites_long <- pivot_longer(clock_rate_2_sites, seq(1,ncol(clock_rate_2_sites)))
clock_rate_2_sites_long$data = "2 sites"

#order columns
clock_rate_all_sites_long <- mutate(clock_rate_all_sites_long,name = fct_relevel(name,ordered_names_all))
clock_rate_2_sites_long <- mutate(clock_rate_2_sites_long, name = fct_relevel(name, ordered_names_2))

clock = rbind(clock_rate_all_sites_long, clock_rate_2_sites_long)

# Test the average of the medians

## Among all sites
tapes_4x = ordered_names_all[2:4]
tapes_5x = ordered_names_all[5:length(ordered_names_all)]

medians_tapes4x_allsites = apply(X = clock_rate_all_sites[, tapes_4x], MARGIN = 2, FUN = median)
mean_tapes4x_allsites = mean(medians_tapes4x_allsites)
mean_tapes4x_allsites
sd(medians_tapes4x_allsites)

medians_tapes5x_allsites = apply(X = clock_rate_all_sites[, tapes_5x], MARGIN = 2, FUN = median)
mean_tapes5x_allsites = mean(medians_tapes5x_allsites)
mean_tapes5x_allsites
sd(medians_tapes5x_allsites)

## Among 2 sites
medians_tapes4x_2sites = apply(X = clock_rate_2_sites[, tapes_4x], MARGIN = 2, FUN = median)
mean_tapes4x_2sites = mean(medians_tapes4x_2sites)
mean_tapes4x_2sites
sd(medians_tapes4x_2sites)

medians_tapes5x_2sites = apply(X = clock_rate_2_sites[, tapes_5x], MARGIN = 2, FUN = median)
mean_tapes5x_2sites = mean(medians_tapes5x_2sites)
mean_tapes5x_2sites
sd(medians_tapes5x_2sites)


# --- 1. Define the color map and legend labels ---
tape_color_map <- c(
  "TTCACGTA" = "#8B0000", 
  "ATGGTAAG" = "#5CA17D", 
  "TTGAGGTG" = "#5CA17D", 
  "TGGTTTTG" = "#FFFF00", 
  "TGCGATTT" = "#5CA17D",
  "TGGACGAC" = "#FFFF00",
  "TTAGATTG" = "#5CA17D",
  "GTAAAGAT" = "#5CA17D",
  "TTTCGTGA" = "#FFFF00",
  "TAGATTTT" = "#5CA17D",
  "ATTTGGTT" = "#5CA17D",
  "GCAGGGTG" = "#5CA17D",
  "ATTTATAT" = "#5CA17D"
)



legend_labels <- c(
  "TTCACGTA" = "2X Tape", 
  "TGGTTTTG" = "4X Tape", 
  "ATGGTAAG" = "5X Tape"
)

p_clock_pos_2 <- ggplot(clock_rate_2_sites_long, aes(x = name, y = value, fill = name)) +
  facet_grid(data ~ .) +
  
  # Background Gridlines (Major Y only)
  # Adding this before geoms so lines are behind the violins
  theme_classic() + 
  
  geom_violin(draw_quantiles = c(0.5), linewidth = 0.2, trim = TRUE) +
  
  scale_fill_manual(
    values = tape_color_map,
    breaks = names(legend_labels),
    labels = legend_labels
  ) +
  
  

  labs(y = expression("Posterior editing rate [" * d^-1 * "]"),x="Tape") +
  
  theme(
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    panel.grid.major.y = element_line(color = "#ebebeb", linewidth = 0.2), 
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1, size = 8),
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 0.5),
    strip.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5),
    panel.spacing = unit(1, "lines")
  )

p_clock_pos_2 

  p_clock_pos_all <- ggplot(clock_rate_all_sites_long, aes(x = name, y = value, fill = name)) +
  facet_grid(data ~ .) +
  
  theme_classic() + 
  
  geom_violin(draw_quantiles = c(0.5), linewidth = 0.2, trim = TRUE) +
  
  scale_fill_manual(
    values = tape_color_map,
    breaks = names(legend_labels),
    labels = legend_labels
  ) +
  
  
  
  labs(y = expression("Posterior editing rate [" * d^-1 * "]"),x="Tape") +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    panel.grid.major.y = element_line(color = "#ebebeb", linewidth = 0.2), 
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1, size = 8),
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 0.5),
    strip.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5),
    panel.spacing = unit(1, "lines")
  )

p_clock_pos_all 

both <- cowplot::plot_grid(p_clock_pos_2,p_clock_pos_all,nrow=2)

ggsave(figure_path, plot = both,height = 14.28, units = "cm", dpi = 300)

