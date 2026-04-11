## ---------------------------
##
## Script name: create_supp_fig_5
##
## Purpose of script: Plot the number of edits introduced at different sites
## based on data in https://doi.org/10.1038/s41586-022-04922-8
##
## Author: Sophie Seidel
##
## Date Created: 2022-08-22
##
## Copyright (c) Sophie Seidel, 2022
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes: data obtained from https://github.com/shendurelab/DNATickerTape/
##
##
## ---------------------------

library(ggplot2)
library(dplyr)
library(reshape2)

figure_path = "plots/supp_fig_5.pdf"
text_size_max = 10
text_size_min = 8
line_thickness = 0.2

# --- Load and Preprocess ---
edit_table_by_5 = readRDS("../figure_3/data/edit_table_filtered.RDS")
write.csv(edit_table_by_5,"supp_fig_5.csv")


# Melt the table so that Site1, Site2, etc. become a single column
# This assumes your site columns are named "Site1", "Site2", "Site3", "Site4", "Site5"
edits_melted <- melt(edit_table_by_5, 
                     id.vars = c("TargetBC", "Cell"), 
                     measure.vars = c("Site1", "Site2", "Site3", "Site4", "Site5"),
                     variable.name = "Sites") %>%
  filter(!is.na(value)) %>% # Equivalent to value != "None" for NA-based tables
  mutate(TargetBC = factor(TargetBC, levels = sort(unique(TargetBC))),
         TargetBC_new = TargetBC,
         SiteNum = as.numeric(gsub("Site", "", Sites)))

# --- Plotting ---
edited_sites_per_tape <- ggplot(edits_melted, aes(x=SiteNum, col=factor(TargetBC_new), group=TargetBC_new)) +
  geom_point(stat = "count", size = 0.5) + 
  geom_line(stat = "count", linewidth = 0.3) +
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

# Display and Save
edited_sites_per_tape
ggsave(plot = edited_sites_per_tape, filename = figure_path, width = 180, height = 100, units = "mm",dpi = 300)


