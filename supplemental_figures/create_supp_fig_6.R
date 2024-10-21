## ---------------------------
##
## Script name: plot_treespace
##
## Purpose of script: run treespace to calculate and plot wRF distances between SciPhy and UPGMA
##
## Author: Antoine Zwaans
##
## Date Created: 2023-11-03
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## -----------------------------


## set output directory for the plot
pic_dir = "plots/"

# set figure settings
text_size=10
label_size = 12


## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(treespace)

##sourcing helper functions
source("plot_trees_in_2d_with_likelihoods.R")

## Load tree data and extract sample nbrs
sciphy_trees <- ape::read.nexus(file = "../figure_4/inference_output/clockPerTarget_sampling_DataSet1_3000000.trees")


sample_nr_tree <- as.numeric(unlist(strsplit(names(sciphy_trees),"_"))[seq(2,2*length(sciphy_trees),by=2)])

## Load log data and extract likelihood values corresponding to the matching sample nrs
log <- read.table("../figure_4/inference_output/clockPerTarget_sampling_DataSet1_3000000.log", header = T)


step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- log$Sample[2] - log$Sample[1]

#resample he log file at the same frequency
subsampled_log <- log[seq(1,length(log$Sample),by=step_tree/step_log),]

#resample the same number
subsampled_log <- subsampled_log[1:min(length(sciphy_trees),length(subsampled_log$Sample)),]

#check that all tree sample nrs and likelihood sample nrs match
which(subsampled_log$Sample != sample_nr_tree)

#extract likelihood values
tree_likelihood <- subsampled_log$likelihood

#if needed check low likelihood sciphy_trees and remove them (potential remnants of burnin)
#plot(tree_likelihood)

#get the upgma tree corresponding to the dataset and relabel the tips to match BEAST tree
upgma <- ape::read.tree(file = "../figure_4/inference_output/UPGMAtree_1000.txt")
cell_ids <- read.csv(header = F, file = "../figure_4/inference_output/UPGMAtree_1000_cell_names_.txt")
cell_ids$numeric_label <- 0:999
cell_ids_sorted <- cell_ids[match(upgma$tip.label, cell_ids$V1), ]
upgma$tip.label <- as.character(cell_ids_sorted$numeric_label)

#create a upgma scaled by the median posterior tree height.
median_posterior_height <- median(log[,"treeHeight.t.alignment"])
upgma_rescaled <- upgma
upgma_height <- tree_height_calc(upgma)
upgma_rescaled$edge.length <- upgma_rescaled$edge.length * (median_posterior_height/upgma_height)

#get the MCC tree 
MCC <- ape::read.nexus(file = "../figure_4/inference_output/MCC_median_heights_clockPerTarget_sampling_DataSet1_3000000.tree")
UPGMA_SCIPHY <- ape::read.nexus(file = "../figure_4/inference_output/combined_1000_UPGMA_medianPosteriorHeight_estimateBranchLengths_infer_rho_sampling_MCC.tree")

#create a list of trees to analyse
all_trees <- sciphy_trees

#append upgma trees and MCC to trees list
all_trees <-c(all_trees, upgma) 
all_trees <-c(all_trees, MCC) 
all_trees <- c(all_trees,upgma_rescaled) 

## ---------------------------------------------------------------------
## Place all all_trees in 2d using the weighted Robinson Foulds (wRF) metric
## ---------------------------------------------------------------------

res_wrf <- treespace(all_trees, nf = 12, method = "wRF")

# Plot wRF scenario
tree_df_wrf <- res_wrf$pco$li
plot1_2 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",1,2) + theme(axis.text = element_blank(),axis.ticks = element_blank(), legend.position ="none") + xlab("Coordinate 1") + ylab("Coordinate 2")
plot1_3 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",1,3) + theme(plot.title = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(), legend.position ="none") + xlab("Coordinate 1") + ylab("Coordinate 3") 
plot1_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",1,4) + theme(plot.title = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(), legend.position ="none") + xlab("Coordinate 1") + ylab("Coordinate 4")
plot2_3 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",2,3) + theme(plot.title = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(), legend.position ="none") + xlab("Coordinate 2") + ylab("Coordinate 3")
plot2_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",2,4) + theme(plot.title = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(), legend.position ="none") + xlab("Coordinate 2") + ylab("Coordinate 4")
plot3_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",3,4) + theme(plot.title = element_blank(), axis.text = element_blank(),axis.ticks = element_blank() ) + xlab("Coordinate 3") + ylab("Coordinate 4")

blank_plot <- ggplot() + theme_void()
combined <- cowplot::plot_grid(plot1_2,plot1_3,plot1_4,blank_plot,plot2_3,plot2_4,blank_plot,blank_plot,plot3_4,nrow = 3,ncol=3)
ggsave(paste0(pic_dir,"supp_fig_6.png"), combined, width = 60, height = 60, units = "cm", dpi = 300)

