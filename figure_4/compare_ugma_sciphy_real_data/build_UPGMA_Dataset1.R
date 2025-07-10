## ---------------------------
##
## Script name:  build_UPGMA_Dataset1
##
## Purpose of script: create upgma from cell culture data subsampled alignments as provided for
## BEAST2.
##
## Author: Antoine Zwaans & Sophie Seidel
##
## Date Created: 2023-06-14
##
## Copyright (c) Antoine Zwaans & Sophie Seidel, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## ---------------------------
##
## Notes: Based on J.Choi's script in DNA typewriter repo.
##
##
## ---------------------------


## load up the packages we will need:
require(tidyverse)
require(data.table)
library(ape)
library(phangorn)
library(tidyverse)
library(ggdendro)
library(dendextend)
library(phytools)

source("upgma_helper_functions.R")

## input sthe filtered set of ccells/sequencies
edit_table_filtered = readRDS("../figure_3/data/edit_table_filtered.RDS")

# sample the same 1000 cells as in Dataset1
set.seed(1)
cell1000 <- sample(unique(edit_table_filtered$Cell),1000)


sample_1000 <- c()
for(i in cell1000) {
  sample_1000 <- rbind(sample_1000,edit_table_filtered[edit_table_filtered$Cell == i,])
}

#create a UPGMA tree for that sample, using the algorithm written by J.Choi
sample_1000 = mark_unedited_sites(sample_1000)
sample_1000 = set_non_contracted_sites_to_na(sample_1000)

# nSites is 13 targetBCs * 5 sites = 65 sites; minus the contracted sites,
# i.e. 1x 3 contracted (2xTape) and 3x 1 contracted (4xTape)
tree = build_upgma_tree(edit_table = sample_1000, nSites = 59)

## we want to scale the tree such that is is smaller than 25 in height.
tree_1000_height <- tree_height_calc(tree)
tree$edge.length <- tree$edge.length * (24.999/tree_1000_height)

## save UPGMA as txt file
ape::write.tree(tree, file='inference_output/UPGMAtree_1000.txt')
write(cell1000,"inference_output/UPGMAtree_1000_cell_names_.txt")

