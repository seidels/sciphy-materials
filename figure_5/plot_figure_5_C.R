## ---------------------------
##
##
## Purpose of script: Visualise edit distribution in mGASv2 dataset
##
## Author: Sophie Seidel
##
## Date Created: 2024-05-20
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


## load up the packages we will need:  

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
## ---------------------------

## load up our functions into memory
source("processing_scripts/useful_scripts_across_categories.R")

## ---------------------------
#figure settings
text_size = 10

## ---------------------------
# Data preprocessing

# input file
filtered_dat_file = "processed_data/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS"
filtered_dat = readRDS(filtered_dat_file)

# sort targetBCs alphabetically and keep record of the ordering
targetbcs_sort = data.frame(tape_names = paste0(1:8), tape_identifier = sort(unique(filtered_dat$TargetBC)))
write.csv(x = targetbcs_sort, file = "tape_seq_identifier_map.csv", quote = F, row.names = F)

# reshape dataframe for plotting
edits_melted = melt(filtered_dat, id.vars = c("TargetBC", "Cell"), variable.name = "Sites")
edits_melted = edits_melted[edits_melted$value != "None", ]
head(edits_melted)

# Add new target bc/ tape names
edits_melted <- edits_melted %>%
  mutate(TargetBC = factor(TargetBC, levels = sort(unique(TargetBC)))) %>%
  mutate(TargetBC_new = paste0(as.numeric(TargetBC)))


edits_melted <- edits_melted %>%
  mutate(SiteNum = as.numeric(gsub("Site", "", Sites)))

## ---------------------------
##### Ask interesting questions ;P

## What is the number of edited sites? (Larger number of edited sites per targetbc --> more info)
head(filtered_dat)

g =
  ggplot(edits_melted, 
         aes(x=SiteNum, col=TargetBC_new, group=TargetBC_new)) +
  geom_point(stat = "count")+
  geom_line(stat = "count")+
  #geom_bar()+
  scale_color_brewer(palette = "Dark2") +
  theme_minimal()+ 
  ylab("Insert count") +
  xlab("Sites")+
  #coord_cartesian(ylim = c(0, 100))+
  theme(text = element_text(size = text_size), 
        legend.position = "right", legend.title = element_blank()) #axis.text.x = element_blank()

g

ggsave( filename = "plots/insert_count_per_site.png", plot = g, width = 6.8, height = 7.2, units = "cm", dpi = 300)

