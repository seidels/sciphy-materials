## Script name: plot_mGASv2_summary_full.R
## Purpose: Final summary figure for mGASv2 analysis
## Dimensions: 180mm x 185mm | Font: 5-7pt

# --- 1. Global Settings & Libraries ---
text_size_max = 7
text_size_min = 5
line_thickness = 0.2

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)
library(reshape2)
library(jsonlite)
library(ggtree)
library(treeio)
library(pammtools)

# --- 2. Data Loading & Processing ---

# A. Clock Rates
typewriter_file <- "inference_output/4-mGASv2-skyline-ou-40K-3-seeds.log"
typewriter <- read.table(typewriter_file, header = T)

clock_rate <- typewriter[, startsWith(colnames(typewriter), "clockRate_")]
order_of_tbcs_in_alignment = c("AGGCTAATTCCC", "TAACGAAGATTT", "AACTGAATGTTT",
                               "GTATAAAGTTTG", "TTGATAACGTGA", "GTTGAAAGGTGA", 
                               "GTACAAAATAGT", "TTCACAAATTTA") 

tape_to_int_map = read.csv(file = "tape_seq_identifier_map.csv")
numeric_names = unname(sapply(order_of_tbcs_in_alignment, function(x){
  which(tape_to_int_map$tape_identifier == x)
}))
names(clock_rate) <- numeric_names
clock_rate <- bind_cols(clock_rate, Prior = rlnorm(nrow(clock_rate), meanlog = -2, sdlog = 0.5))
clock_rate_long <- pivot_longer(clock_rate, everything())

# B. Site Editing
filtered_dat <- readRDS("processed_data/mGASv2_Lane2_CellByTape_filtered_for_8barcodes_OG.RDS")
edits_melted <- melt(filtered_dat, id.vars = c("TargetBC", "Cell"), variable.name = "Sites") %>%
  filter(value != "None") %>%
  mutate(TargetBC = factor(TargetBC, levels = sort(unique(TargetBC))),
         TargetBC_new = as.numeric(TargetBC),
         SiteNum = as.numeric(gsub("Site", "", Sites)))

# C. Growth Rates (SciPhy + Merle et al.)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.", 1:3)] - typewriter_mcmc[,paste0("deathRate.", 1:3)]
HPD <- HPDinterval(growth)
growth_hpd <- data.frame(median = as.numeric(apply(growth, 2, median)),
                         hpd_low = as.numeric(HPD[,"lower"]),
                         hpd_up = as.numeric(HPD[,"upper"]))
growth_combined <- rbind(growth_hpd, growth_hpd[3,]) 
growth_combined$t <- c(0.0, 4, 7.5, 11)
growth_combined$tree <- "SciPhy"

# Merle et al. JSON processing
merle_raw <- read_json(path="data/MerleEtAl_Extended_Data_Fig_2b_raw_data.json", simplifyVector = TRUE)$data
rate <- c()
for(i in 0:22) { 
  for(j in c(2,5)) {
    if(j==2) rate <- c(rate, log(merle_raw$N[i*5 + j]/merle_raw$N[i*5 + j-1]))
    if(j==5) rate <- c(rate, (1/3)*log(merle_raw$N[i*5 + j]/merle_raw$N[i*5 + j-3]))
  }
}
merle_time <- data.frame(T = rep(c(6.0, 7.0), 23), growth = rate)
m_mean <- aggregate(growth ~ T, merle_time, mean)$growth
m_sd <- aggregate(growth ~ T, merle_time, sd)$growth
merle_et_al_bins <- data.frame(t = c(6.0, 7.0, 11), 
                               mean = c(m_mean, m_mean[2]), 
                               low = c(m_mean - m_sd, (m_mean - m_sd)[2]), 
                               high = c(m_mean + m_sd, (m_mean + m_sd)[2]), 
                               tree = "Merle et al., 2024")

# D. Tree
tree <- read.beast("inference_output/4-mGASv2-skyline-ou-40K.1000Kresampled-3-seeds.CCD0.defaultHeights.txt")
origin_length <- 11 - round(max(tree@data$height, na.rm=T), 3)

# --- 3. Individual Plots ---
dark2_cols <- scales::brewer_pal(palette = "Dark2")(8)
names(dark2_cols) <- as.character(sort(numeric_names))
custom_palette <- c(dark2_cols, "Prior" = "#E1E1F7")

p_clock <- ggplot(clock_rate_long, aes(x=name, y=value, fill=name)) +
  geom_violin(draw_quantiles = 0.5, linewidth = line_thickness, adjust = 1.8, trim = FALSE) +
  scale_fill_manual(values = custom_palette) + 
  theme_classic(base_size = text_size_min) +
  labs(y = expression("Editing rate [" * d^-1 * "]"), x = "Estimates per tape") +
  theme(legend.position = "none", 
        axis.title = element_text(size = text_size_max),
        axis.text = element_text(size = text_size_min),
        axis.line = element_line(linewidth = line_thickness), 
        axis.ticks = element_line(linewidth = line_thickness))

p_editing <- ggplot(edits_melted, aes(x=SiteNum, col=factor(TargetBC_new), group=TargetBC_new)) +
  geom_point(stat = "count", size = 0.5) + geom_line(stat = "count", linewidth = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  theme_classic(base_size = text_size_min) +
  labs(y = "Insert count", x = "Sites", color = "Tape") +
  theme(legend.key.spacing.y = unit(0.02, "cm"),
        legend.title = element_text(size = text_size_min),
        legend.text = element_text(size = text_size_min),
        legend.position = "right", 
        axis.title = element_text(size = text_size_max),
        axis.text = element_text(size = text_size_min),
        axis.line = element_line(linewidth = line_thickness), 
        axis.ticks = element_line(linewidth = line_thickness))

p_growth <- ggplot(growth_combined) +  
  geom_step(aes(x=t, y=median, col=tree), linewidth = 0.5) + 
  geom_stepribbon(aes(x=t, ymin = hpd_low, ymax=hpd_up, fill=tree), alpha = 0.1) +
  geom_step(data=merle_et_al_bins, aes(x=t, y=mean, col=tree), linewidth = 0.5, linetype="dashed") +
  geom_stepribbon(data=merle_et_al_bins, aes(x=t, ymin=low, ymax=high, fill=tree), alpha = 0.1) +    
  scale_color_manual(values = c("SciPhy"="red", "Merle et al., 2024"="black")) +
  scale_fill_manual(values = c("SciPhy"="red", "Merle et al., 2024"="black")) +
  theme_classic(base_size = text_size_min) + 
  labs(y = expression("Growth rate [" * d^-1 * "]"), x = "Time [d]") +
  theme(legend.position = c(0.8, 0.8), 
        legend.title = element_blank(), 
        legend.text = element_text(size = text_size_min),
        axis.title = element_text(size = text_size_max),
        axis.text = element_text(size = text_size_min),
        axis.line = element_line(linewidth = line_thickness), 
        axis.ticks = element_line(linewidth = line_thickness)) + 
  scale_x_continuous(breaks=c(0,4,6,7,7.5,10,11)) 

p_tree <- ggtree(tree, root.position = origin_length, size=0.15, color="darkgrey") + 
  theme_tree2() + geom_rootedge(rootedge = origin_length, size=0.15) + 
  geom_nodepoint(aes(size=posterior)) + scale_size(range = c(0.1, 2)) +
  vexpand(.05, 1) +  vexpand(.05, -1) +
  scale_x_continuous(breaks = c(0, 4, 6, 7, 7.5, 11)) + 
  labs(x = "Time [d]", size="Posterior support") +
  theme(legend.position = "top",
        legend.text = element_text(size = text_size_min),
        legend.title = element_text(size = text_size_min),
        text = element_text(size = text_size_min), 
        axis.text.x = element_text(size = text_size_min),
        axis.title.x = element_text(size = text_size_max),
        axis.line.x = element_line(linewidth = line_thickness),
        axis.ticks.x = element_line(linewidth = line_thickness))

# --- 4. Assembly ---
col_bc <- plot_grid(p_clock, p_editing, ncol = 2, labels = c("B", "C"), label_size = text_size_max)
final_plot <- plot_grid(ggdraw(), col_bc, p_growth, p_tree, ncol = 1, 
                        labels = c("A", "", "D", "E"), label_size = text_size_max, rel_heights = c(1, 1, 1, 2.0))

ggsave("plots/figure_mGASv2_GUIDELINE_ADJUSTED_correct_font_size.pdf", final_plot, width = 180, height = 185, units = "mm", dpi = 300, device = cairo_pdf)
