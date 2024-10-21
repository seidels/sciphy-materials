## ---------------------------
##
## Script name: plot_trees_in_2d_with_likelihoods
##
## Purpose of script: Show similarity of UPGMA tree to
## the trees in the posterior set with respect to tree
## topology only (measured by Robinson Foulds) and with
## respect to tree topology and branch lengths (weighted
## Robinson Foulds.
##
## Based on Sophie Seidel's plot_trees_in_2d.R script
##
## Date Created: 2023-07-21, Antoine Zwaans
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##

## Load required packages
library(treespace)
library(adegraphics)
library(ggplot2)
library(rgl)
library(phangorn)

## ---------------------------

## Define functions

# Function to plot tree distances colored by likelihood
plot_tree_distances_topology_metric <- function(tree_df,tree_likelihood,metric_name,dim1,dim2) {
 tree_df <-  tree_df
 plot_title = paste(metric_name,"metric space")
  dim_1 <- paste0("A",dim1)
  dim_2 <- paste0("A",dim2)
   
  plot <- ggplot(tree_df[1: (nrow(tree_df) - 2),], aes_string(x = dim_1, y = dim_2)) +
    geom_point(aes(col = tree_likelihood), size = 2, alpha = 0.5) + 
    scale_color_gradient(low = "#E1E1F7", high = "#060647",name = "Likelihood") + 
    xlab(dim_1) + ylab(dim_2) + theme_bw(base_family = "")+ 
    theme(legend.position = c(0.6,1.43), text = element_text(size = 10)) +
    geom_point(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),colour="red",size = 2)+ 
    geom_text(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),label="UPGMA",hjust = -0.2,vjust =1,size=2.5)+ 
    geom_point(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),colour="black",size = 2)+ 
    geom_text(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),label="MCC",hjust = -0.3,vjust =1,size=2.5) +
    ggtitle(plot_title)
  
  return(plot)
}

#red color grdient, replace with
#scale_color_gradient(low = "#E0C7C7", high = "#990000",name = "Likelihood")

#green gradient
#scale_color_gradient(low = "#A3D6AA", high = "#075713",name = "Likelihood")
#for the figure (dot size):
#size=2


plot_tree_distances_topology_metric_MAP <- function(tree_df,tree_likelihood,metric_name,dim1,dim2,index_MAP_tree) {
  tree_df <-  tree_df
  plot_title = paste(metric_name,"metric space")
  dim_1 <- paste0("A",dim1)
  dim_2 <- paste0("A",dim2)
  
  plot <- ggplot(tree_df[1: (nrow(tree_df) - 1),], aes_string(x = dim_1, y = dim_2)) +
    geom_point(aes(col = tree_likelihood), size = 6, alpha = 0.5) + 
    scale_color_gradient(low = "#E1E1F7", high = "#060647",name = "Likelihood")  +  
    xlab(dim_1) + ylab(dim_2) + theme_bw(base_family = "")+ 
    theme(legend.position = c(0.8,0.2), text = element_text(size = 14)) +
    geom_point(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),colour="red",size=5)+ 
    geom_text(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),label="UPGMA",hjust = -0.2,vjust =1,size=6)+ 
    geom_point(data=tree_df[(index_MAP_tree):(index_MAP_tree),], aes_string(x=dim_1,y=dim_2),colour="red",size=5)+ 
    geom_text(data=tree_df[(index_MAP_tree):(index_MAP_tree),], aes_string(x=dim_1,y=dim_2),label="MAP",hjust =1,vjust =-1,size=6) +
    ggtitle(plot_title)
  
  return(plot)
}

plot_tree_distances_topology_metric_MSCC <- function(tree_df,tree_likelihood,metric_name,dim1,dim2) {
  tree_df <-  tree_df
  plot_title = paste(metric_name,"metric space")
  dim_1 <- paste0("A",dim1)
  dim_2 <- paste0("A",dim2)
  
  plot <- ggplot(tree_df[1: (nrow(tree_df) - 3),], aes_string(x = dim_1, y = dim_2)) +
    geom_point(aes(col = tree_likelihood), size = 6, alpha = 0.5) + 
    scale_color_gradient(low = "#E1E1F7", high = "#060647",name = "Likelihood")  +  
    xlab(dim_1) + ylab(dim_2) + theme_bw(base_family = "")+ 
    theme(legend.position = c(0.8,0.2), text = element_text(size = 14)) +
    geom_point(data=tree_df[(nrow(tree_df) - 2):(nrow(tree_df) - 2),], aes_string(x=dim_1,y=dim_2),colour="red")+ 
    geom_text(data=tree_df[(nrow(tree_df) - 2):(nrow(tree_df) - 2),], aes_string(x=dim_1,y=dim_2),label="UPGMA",hjust = 1,vjust =1)+ 
    geom_point(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),colour="red")+ 
    geom_text(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),label="MCC",hjust = 1,vjust =1)+ 
    geom_point(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),colour="red")+ 
    geom_text(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),label="MSCC",hjust =1,vjust =-1) +
    ggtitle(plot_title)
  
  return(plot)
}

plot_tree_distances_branch_lengths_metric <- function(tree_df,tree_likelihood,metric_name,dim1,dim2) {
  tree_df <-  tree_df
  plot_title = paste("Weigthed Robinson Foulds","metric space")
  dim_1 <- paste0("A",dim1)
  dim_2 <- paste0("A",dim2)
  
  plot <- ggplot(tree_df[1: (length(tree_df$A1) - 3),], aes_string(x = dim_1, y = dim_2)) +
    geom_point(aes(col = tree_likelihood), size = 2, alpha = 0.5) + 
    scale_color_gradient(low = "#E1E1F7", high = "#060647",name = "Likelihood")  +  
    xlab("") + ylab("") + theme_bw(base_family = "")+ 
    theme(legend.position = c(0.8,0.2), text = element_text(size = 10)) +
    geom_point(data=tree_df[(nrow(tree_df) - 2):(nrow(tree_df) - 2),], aes_string(x=dim_1,y=dim_2),colour="red",size = 2)+ 
    geom_text(data=tree_df[(nrow(tree_df) - 2):(nrow(tree_df) - 2),], aes_string(x=dim_1,y=dim_2),label="UPGMA",hjust = -0.2,vjust =1)+ 
    geom_point(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),colour="red",size = 2.5)+ 
    geom_text(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),label="UPGMA root",hjust = -0.3,vjust =-0.2) +
    geom_point(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),colour="black",size = 2.5)+ 
    geom_text(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),label="MCC",hjust = -0.2,vjust =1) +
    ggtitle(plot_title)
    
  
  return(plot)
}

tree_height_calc <- function(tree) { 
  start_edge <- 1
  sum_path <- 0
  while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
    sum_path <- sum_path + tree$edge.length[tree$edge[,2] == start_edge]
    start_edge = tree$edge[tree$edge[,2] == start_edge,1]
    
  }
  return(sum_path)
  
}

