## ---------------------------
##
## Script name: create_supp_fig_12
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
library(phytools)

##sourcing helper functions
plot_tree_distances_branch_lengths_metric <- function(tree_df,tree_likelihood,metric_name,dim1,dim2) {
  tree_df <-  tree_df
  plot_title = paste("Weigthed Robinson Foulds","metric space")
  dim_1 <- paste0("A",dim1)
  dim_2 <- paste0("A",dim2)
  
  plot <- ggplot(tree_df[1: (length(tree_df$A1) - 3),], aes_string(x = dim_1, y = dim_2)) +
    geom_point(aes(col = tree_likelihood), size = 0.4, alpha = 0.5) + 
    scale_color_gradient(low = "white", high = "#05290B",name = "Log-likelihood")  +  
    xlab("") + ylab("") + theme_bw(base_family = "")+ 
    theme(legend.position = c(0.8,0.2), text = element_text(size = 10))+
    geom_point(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),colour="orange",size = 2.5)+ 
    geom_text(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),label="UPGMA root",hjust = 0.3,vjust =-0.2,size=2.5) +
    geom_point(data=tree_df[(nrow(tree_df) - 2):(nrow(tree_df) - 2),], aes_string(x=dim_1,y=dim_2),colour="red",size = 2.5)+ 
    geom_text(data=tree_df[(nrow(tree_df) - 2):(nrow(tree_df) - 2),], aes_string(x=dim_1,y=dim_2),label="UPGMA",hjust = 0.2,vjust =1,size=2.5)+ 
    geom_point(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),colour="black",size = 2.5)+ 
    geom_text(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),label="CCD",hjust = -0.2,vjust =1,size=2.5) +
    ggtitle(plot_title) + theme(axis.text = element_text(size=6))
  
  
  return(plot)
}

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
upgma_height <- max(nodeHeights(upgma))
upgma_rescaled$edge.length <- upgma_rescaled$edge.length * (median_posterior_height/upgma_height)

#get the MCC tree 
MCC <- ape::read.nexus(file = "../figure_4/inference_output/CCD_clockPerTarget_sampling_DataSet1_3000000.tree")
#UPGMA_SCIPHY <- ape::read.nexus(file = "../figure_4/inference_output/CCD_UPGMA_SciPhy.tree")

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
write.csv(res_wrf$pco$li,"supp_fig_12.csv")

plot1_2 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weighted Robinson Foulds",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weighted Robinson Foulds",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weighted Robinson Foulds",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_3 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weighted Robinson Foulds",2,3) + xlab("Coordinate 2") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weighted Robinson Foulds",2,4) + xlab("Coordinate 2") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot3_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weighted Robinson Foulds",3,4) + xlab("Coordinate 3") + ylab("Coordinate 4") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

blank_plot <- ggplot() + theme_void()
legend <- get_legend(plot3_4) 
combined <- cowplot::plot_grid(plot1_2+ theme(plot.title = element_blank()),plot1_3 + theme(plot.title = element_blank()),plot1_4+ theme(plot.title = element_blank()),blank_plot,plot2_3 + theme(plot.title = element_blank()),plot2_4 + theme(plot.title = element_blank()),legend,blank_plot,plot3_4+ theme(plot.title = element_blank(),legend.position = "None"),nrow = 3,ncol=3)

#plot with titles
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Weighted Robinson Foulds metric space",
    fontface = 'bold',size = 10,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
combined_with_title <- cowplot::plot_grid(
  title, combined,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave(paste0(pic_dir,"supp_fig_12.pdf"),combined_with_title, width = 14.28, height = 14.28, units = "cm", dpi = 300)

