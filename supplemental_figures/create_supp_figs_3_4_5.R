## ---------------------------
##
## Script name: create_supp_fig_8_9_10_11
##
## Purpose of script: Plot outputs of 2D mapping obtained with TreeDist and assess quality.
## Relies on outputs from run_TreeDist_XX.R, unmapped distances and mapped distance run on cluster.
##
## Author: Antoine Zwaans
##
## Date Created: 2023-09-11
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## -----------------------------

## set directory

setwd("/Users/azwaans/typewriter_analysis/paper_figures/figure_4/")

## set output directory for the plot
pic_dir = "plots/"

# set figure settings
text_size=10
label_size = 12


## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(data.table)
library(ape)
library(TreeDist)
library(gdata)

#########################################################
## Getting likelihood values matching the SciPhy trees ##
#########################################################

##sourcing helper functions
plot_tree_distances_topology_metric <- function(tree_df,tree_likelihood,metric_name,dim1,dim2) {
  tree_df <-  tree_df
  plot_title = paste(metric_name,"metric space")
  dim_1 <- paste0("A",dim1)
  dim_2 <- paste0("A",dim2)
  
  plot <- ggplot(tree_df[1: (nrow(tree_df) - 2),], aes_string(x = dim_1, y = dim_2)) +
    geom_point(aes(col = tree_likelihood), size = 0.4, alpha = 0.5) + 
    scale_color_gradient(low = "white", high = "#05290B",name = "Log-likelihood") + 
    xlab(dim_1) + ylab(dim_2) + theme_bw(base_family = "")+ 
    theme(legend.position = c(0.6,1.43), text = element_text(size = 10)) +
    geom_point(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),colour="red",size = 2)+ 
    geom_text(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),label="UPGMA",hjust = -0.2,vjust =1,size=2.5)+ 
    geom_point(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),colour="black",size = 2)+ 
    geom_text(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),label="CCD",hjust = -0.3,vjust =1,size=2.5) +
    ggtitle(plot_title) + theme(axis.text = element_text(size=6))
  
  return(plot)
}

#### preprocess log to extract the correct likelihood values to plot ####
sciphy_trees <- ape::read.nexus(file = "../figure_4/inference_output/clockPerTarget_sampling_DataSet1_3000000.trees")

sample_nr_tree <- as.numeric(unlist(strsplit(names(sciphy_trees),"_"))[seq(2,2*length(sciphy_trees),by=2)])

## Load log data and extract likelihood values corresponding to the matching sample nrs
log <- read.table("../figure_4/inference_output/clockPerTarget_sampling_DataSet1_3000000.log", header = T)

step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- log$Sample[2] - log$Sample[1]

#resample the log file at the same frequency
subsampled_log <- log[seq(1,length(log$Sample),by=step_tree/step_log),]

#resample the same number
subsampled_log <- subsampled_log[1:min(length(sciphy_trees),length(subsampled_log$Sample)),]

#check that all tree sample nrs and likelihood sample nrs match
which(subsampled_log$Sample != sample_nr_tree)

#extract likelihood values
tree_likelihood <- subsampled_log$likelihood

#############################
## Plotting in RF 2D space ##
############################

#input RF 2D mapping 
RF_treeDist <- read.csv( file = "../figure_4/inference_output/mapping_df_RF_CCD.csv")

#input original distances
distances <- read.csv( file = "../figure_4/inference_output/distances_RF_CCD.csv")

#reformat into lowertriangular matrix
distances <- unlist(distances)
mat <- matrix(NA, ncol=length(sciphy_trees) +2, nrow=length(sciphy_trees) +2) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances
dist <- as.dist(mat)
#calculate trustworthyness x continuity (mapping quality)
txc <- vapply(seq_len(ncol(RF_treeDist)), function(k) {
  newDist <- dist(RF_treeDist[, seq_len(k)])
  MappingQuality(dist, newDist, 10)["TxC"]
}, 0)
plot(txc, xlab = "Dimension")

#plot the mapping quality
png(paste0(pic_dir,"supp_fig_11a.png"),width = 14.28, height = 14.28, units = "cm", res = 300)
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

tree_df_rf <- data.frame(RF_treeDist)
plot1_2_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_3_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",2,3) + xlab("Coordinate 2") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_4_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",2,4) + xlab("Coordinate 2") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot3_4_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",3,4) + xlab("Coordinate 3") + ylab("Coordinate 4") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

blank_plot <- ggplot() + theme_void()
legend <- get_legend(plot3_4_RF) 
combined <- cowplot::plot_grid(plot1_2_RF+ theme(plot.title = element_blank()),plot1_3_RF + theme(plot.title = element_blank()),plot1_4_RF+ theme(plot.title = element_blank()),blank_plot,plot2_3_RF + theme(plot.title = element_blank()),plot2_4_RF + theme(plot.title = element_blank()),legend,blank_plot,plot3_4_RF+ theme(plot.title = element_blank(),legend.position = "None"),nrow = 3,ncol=3)

#plot with titles
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Robinson Foulds metric space",
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

ggsave(paste0(pic_dir,"supp_fig_8.pdf"), combined_with_title, width = 14.28, height = 14.28, units = "cm", dpi = 300)

###########################################
## Plotting in 2D Clustering Information ##
###########################################

CI_treeDist <- read.csv( file = "../figure_4/inference_output/mapping_df_CI_CCD.csv")

#input original distances
distances <- read.csv( file = "../figure_4/inference_output/distances_CI_CCD.csv")

#reformat into lowertriangular matrix
distances <- unlist(distances)
mat <- matrix(NA, ncol=length(sciphy_trees) +2, nrow=length(sciphy_trees) +2) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances
dist <- as.dist(mat)

#calculate trustworthyness x continuity (mapping quality)
txc <- vapply(seq_len(ncol(CI_treeDist)), function(k) {
  newDist <- dist(CI_treeDist[, seq_len(k)])
  MappingQuality(dist, newDist, 10)["TxC"]
}, 0)

#plot the mapping quality
png(paste0(pic_dir,"supp_fig_11b.png"),width = 14.28, height = 14.28, units = "cm", res = 300)
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

#remove the legend, the other plit will have it

tree_df_rf <- CI_treeDist
plot1_2 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",2,3) + xlab("Coordinate 2") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",2,4) + xlab("Coordinate 2") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot3_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",3,4) + xlab("Coordinate 3") + ylab("Coordinate 4") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

blank_plot <- ggplot() + theme_void()
legend <- get_legend(plot3_4) 
combined <- cowplot::plot_grid(plot1_2+ theme(plot.title = element_blank()),plot1_3 + theme(plot.title = element_blank()),plot1_4+ theme(plot.title = element_blank()),blank_plot,plot2_3 + theme(plot.title = element_blank()),plot2_4 + theme(plot.title = element_blank()),legend,blank_plot,plot3_4+ theme(plot.title = element_blank(),legend.position = "None"),nrow = 3,ncol=3)

#plot with titles
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Clustering Information metric space",
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

ggsave(paste0(pic_dir,"supp_fig_9.pdf"),combined_with_title, width = 14.28, height = 14.28, units = "cm", dpi = 300)


#############################################
## Plotting in 2D Phylogenetic Information ##
#############################################

PI_treeDist <- read.csv( file = "../figure_4/inference_output/mapping_df_PI_CCD.csv")

#input original distances
distances <- read.csv( file = "../figure_4/inference_output/distances_PI_CCD.csv")

#reformat into lowertriangular matrix
distances <- unlist(distances)
mat <- matrix(NA, ncol=length(sciphy_trees) +2, nrow=length(sciphy_trees) +2) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances
dist <- as.dist(mat)

#calculate trustworthyness x continuity (mapping quality)
txc <- vapply(seq_len(ncol(PI_treeDist)), function(k) {
  newDist <- dist(PI_treeDist[, seq_len(k)])
  MappingQuality(dist, newDist, 10)["TxC"]
}, 0)

#plot the mapping quality
png(paste0(pic_dir,"supp_fig_11c.png"),width = 14.28, height = 14.28, units = "cm", res = 300)
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

#remove the legend, the other plit will have it

tree_df_rf <- PI_treeDist
plot1_2 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot2_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",2,3) + xlab("Coordinate 2") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",2,4) + xlab("Coordinate 2") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot3_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",3,4) + xlab("Coordinate 3") + ylab("Coordinate 4") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

blank_plot <- ggplot() + theme_void()
legend <- get_legend(plot3_4) 
combined <- cowplot::plot_grid(plot1_2+ theme(plot.title = element_blank()),plot1_3 + theme(plot.title = element_blank()),plot1_4+ theme(plot.title = element_blank()),blank_plot,plot2_3 + theme(plot.title = element_blank()),plot2_4 + theme(plot.title = element_blank()),legend,blank_plot,plot3_4+ theme(plot.title = element_blank(),legend.position = "None"),nrow = 3,ncol=3)

#plot with titles
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Phylogenetic Information metric space",
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

ggsave(paste0(pic_dir,"supp_fig_10.pdf"),combined_with_title, width = 14.28, height = 14.28, units = "cm", dpi = 300)

