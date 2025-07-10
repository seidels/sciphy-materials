## ---------------------------
##
## Script name: TreeDist 
##
## Purpose of script: Running tree dist distance and generate a dataframe for plotting tree posterior in tree space
##
## Author: Antoine Zwaans
##
## Date Created: 2023-09-11
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

library(tidyverse)
library(data.table)
library(ape)
library(TreeDist)

#load sciphy trees
sciphy_trees <- ape::read.nexus(file = "clockPerTarget_sampling_DataSet1_3000000.trees")
all_trees <- sciphy_trees

#load and add upgma tree
upgma <- ape::read.tree(file = "UPGMAtree_1000.txt")
cell_ids <- read.csv(header = F, file = "UPGMAtree_1000_cell_names_.txt")
cell_ids$numeric_label <- 0:999
cell_ids_sorted <- cell_ids[match(upgma$tip.label, cell_ids$V1), ]
upgma$tip.label <- as.character(cell_ids_sorted$numeric_label)
all_trees <- c(all_trees,upgma)

#get the MCC tree and add to the all_trees list
MCC <- ape::read.nexus(file = "MCC_median_heights_clockPerTarget_sampling_DataSet1_3000000.tree")
all_trees <- c(all_trees,MCC)

distances <- RobinsonFoulds(all_trees)
distances_df <- data.frame(distances)
write_csv(distances_df,"distances_RF.csv")
mapping <- cmdscale(distances, k = 12)
colnames(mapping) <- paste0(rep("A",12),1:12)
mapping_df <- data.frame(mapping)
write_csv(mapping_df,"mapping_df_RF.csv")
