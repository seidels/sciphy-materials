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
source("../figure_5/processing_scripts/useful_scripts_across_categories.R")

## ---------------------------
#figure settings
text_size = 10

## ---------------------------
# Data preprocessing

# input file
filtered_dat_file = "./../figure_5/data/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS"
filtered_dat = readRDS(filtered_dat_file)

# sort targetBCs alphabetically and keep record of the ordering
targetbcs_sort = data.frame(tape_names = paste0(1:8), tape_identifier = sort(unique(filtered_dat$TargetBC)))

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
## What is the number of unique sequences? (Larger number of unique sequences --> more info)


# Reduce to cells with different sequences
filtered_dat_unique_cells = unique(filtered_dat[, 2:7])

nr_unique_sequences_per_targetbc = sapply(targetbcs_sort$tape_identifier, function(x){
  nrow(filtered_dat_unique_cells[which(filtered_dat_unique_cells$TargetBC ==x), ])
})

unique_seqs_per_targetbc = data.frame(targetbc = targetbcs_sort$tape_identifier, 
                                      targetbc_new = targetbcs_sort$tape_names,
                                      unique_seqs = nr_unique_sequences_per_targetbc)

g = ggplot(unique_seqs_per_targetbc, aes(x=targetbc_new, y=unique_seqs, fill=targetbc_new)) + 
  geom_col()+
  ylab("Unique seq count")+
  xlab("Tapes")+
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal()+
  theme(legend.position = "none", 
        axis.title.y = element_text(hjust = 0.2))
g
ggsave( filename = "plots/supp_fig_13.pdf", plot = g, width = 7.5, height = 5, units = "cm", dpi = 300)
