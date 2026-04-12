## ---------------------------
##
## Script name: Robustness to error/dropout
##
## Purpose of script: Showcase robustness of SciPhy to dropout silencing
##
## Author: Antoine Zwaans
##
## Date Created: 2025-10-01
##
## ---------------------------

#functions to plot half violin plots side to side
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



require(tidyverse)
require(data.table)
library(ape)
library(TreeSim)
library(phylobase)
library(phangorn)
library(cowplot)
text_size=10

#source("useful_scripts_across_categories.R")

create_distances_to_truth <- function(dataset_name) {
sampling <- dataset_name
setwd(paste0("supp_fig_lossy_data/final_validations_manuscript_",dataset_name))

distance_SCIPHY_truth_PI_complete <- c()
distance_SCIPHY_truth_PI_filtered <- c()

distance_SCIPHY_truth_wrf_complete <- c()
distance_SCIPHY_truth_wrf_filtered <- c()

distance_SCIPHY_truth_ks_filtered <- c()
distance_SCIPHY_truth_ks_complete <- c()

#collect ground truth/simulation parameters
clock_rates <- c()
loss_rates <- c()
loss_probs <- c()
filtered_tree_heights <- c()
tree_heights <- c()
filtered_ntips <- c()
unfiltered_ntips <- c()

for(SEED in 1:100) {
  
  #parse the trees
  print(paste0("seed: ",SEED))

  #check if file is there 
  if(file.exists(paste0("supp_fig_lossy_data/ccd_tree_without_error.",SEED,".txt"))) {
  
  #read-in both trees
  MCC_tree <- ape::read.nexus(paste0("supp_fig_lossy_data/ccd_tree.",SEED,".txt"))
  MCC_tree_filtered <- ape::read.nexus(paste0("supp_fig_lossy_data/ccd_tree_without_error.",SEED,".txt"))
  
  #only read in if tree is at least a cherry, read in simulation parameters
  if(MCC_tree_filtered$Nnode > 1) {
    
    file_name <-  paste0("supp_fig_lossy_data/parameters/simParams_",SEED,".csv")
    parameters <- read.csv(file_name)
    clock_rates <- c(clock_rates,parameters$x[14])
    loss_rates <- c(loss_rates,parameters$x[15] * parameters$x[14])
    loss_probs <- c(loss_probs,parameters$x[16])
 
  #parse the true tree
  file_name_true_tree <- paste0("supp_fig_lossy_data/data/simulate_with_both_sampling=",dataset_name,".",SEED,".newick")
  true_tree <- readLines(file_name_true_tree)
  true_tree <- str_remove_all(true_tree,"\\[.........\\]")
  true_tree <- str_remove_all(true_tree,"\\[..........\\]")
  true_tree <- str_remove_all(true_tree,"\\[...........\\]")
  true_tree <- paste0(true_tree,";")
  true_tree <- ape::read.tree(text = true_tree)
  
  #filter the tree to construct the ground truth tree for the filtered data
  true_tree_filtered <- true_tree
  MCC_tree_filtered$tip.label <- sort(true_tree$tip.label)[as.numeric(MCC_tree_filtered$tip.label)+1]
  
  for(true_tree_tip in true_tree_filtered$tip.label) { 
    
    if(!(true_tree_tip %in% MCC_tree_filtered$tip.label)) {
      print("trying to remove tip ")
      true_tree_filtered <- ape::drop.tip(true_tree_filtered,true_tree_tip)
      
    }
    
    
  }
  
  distance_SCIPHY_truth_PI_filtered <- c(distance_SCIPHY_truth_PI_filtered,TreeDist::PhylogeneticInfoDistance(MCC_tree_filtered,true_tree_filtered,normalize = TRUE))
  distance_SCIPHY_truth_PI_complete <- c(distance_SCIPHY_truth_PI_complete,TreeDist::PhylogeneticInfoDistance(true_tree, MCC_tree,normalize = TRUE))

  distance_SCIPHY_truth_wrf_filtered <- c(distance_SCIPHY_truth_wrf_filtered,wRF.dist(true_tree_filtered,MCC_tree_filtered,normalize = TRUE))
  distance_SCIPHY_truth_wrf_complete <- c(distance_SCIPHY_truth_wrf_complete,wRF.dist(true_tree,MCC_tree,normalize = TRUE))
  
  distance_SCIPHY_truth_ks_filtered <- c(distance_SCIPHY_truth_ks_filtered,ks.test(true_tree_filtered$edge.length,MCC_tree_filtered$edge.length)$statistic)
  distance_SCIPHY_truth_ks_complete <- c(distance_SCIPHY_truth_ks_complete,ks.test(true_tree$edge.length,MCC_tree$edge.length)$statistic)
  

  }
  }
  
}

return(list(complete_PI=distance_SCIPHY_truth_PI_complete,filtering_PI=distance_SCIPHY_truth_PI_filtered,complete_wRF=distance_SCIPHY_truth_wrf_complete,filtering_wRF=distance_SCIPHY_truth_wrf_filtered,heritable=loss_rates,drop=loss_probs,complete_ks=distance_SCIPHY_truth_ks_complete,filtering_ks=distance_SCIPHY_truth_ks_filtered))
}


### 0.00003
distances_0.00003 <- create_distances_to_truth("0.00003")
distances<- c(distances_0.00003$complete_PI,distances_0.00003$filtering_PI,distances_0.00003$complete_wRF,distances_0.00003$filtering_wRF)
metrics<- c(rep("PI distances",length(c(distances_0.00003$complete_PI,distances_0.00003$filtering_PI))),rep("wRF distances",length(c(distances_0.00003$complete_wRF,distances_0.00003$filtering_wRF))))
modellings<- c(rep("Complete data",length(distances_0.00003$complete_PI)),rep("Filtered lossy data",length(distances_0.00003$filtering_PI)),rep("Complete data",length(distances_0.00003$complete_wRF)),rep("Filtered lossy data",length(distances_0.00003$filtering_wRF)))
dataframe_0.00003 <- data.frame(distance = distances,metric=metrics,sampling=rep("rho=0.00003",length(distances)),data_regime=modellings)

###
# combined loss parameter
###

g_distances = ggplot(dataframe_0.00003, aes( x=factor(metric), y=distance, fill=data_regime)) +
  geom_split_violin(width=1.2,adjust=1.3) +
  geom_boxplot(width=0.1,  alpha=0.2, outlier.shape = NA,coef = 0) +
  xlab("")+
  ylab("Normalised distance to true tree")+
  theme_classic()+ 
  theme(strip.text.x = element_text(size = 11),plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = c(0.7,0.7))


### plot relationship between loss probability and distances
distances<- c(distances_0.00003$filtering_PI,distances_0.00003$filtering_wRF)
heritable_rates<- c(distances_0.00003$heritable,distances_0.00003$heritable)
dropout_rates <- c(distances_0.00003$drop,distances_0.00003$drop)
metrics<- c(rep("PI distances",length(distances_0.00003$filtering_PI)),rep("wRF distances",length(distances_0.00003$filtering_wRF)))

dataframe_0.00003 <- data.frame(distance = distances,metric=metrics,sampling=rep("rho=0.00003",length(distances)),heritab=heritable_rates,droping=dropout_rates)


g_distances_scatter_overall = ggplot(dataframe_0.00003, aes( x=exp(-25*heritab)*droping +  1 - exp(-25*heritab) , y=distance)) +
  geom_point() +
  xlab("Combined tape loss probability")+
  ylab("Normalised distance to true tree")+
  theme_classic()+ 
  facet_wrap(facets=~metric) + 
  theme(strip.text.x = element_text(size = 11),plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = c(0.5,0.5))

both_plots_1 <- plot_grid(g_distances,g_distances_scatter_overall,ncol=1)


ggsave(plot = both_plots_1, "plots/supp_fig_23.pdf",
       width = 14.28, height = 14.28, units = "cm")


###
# Individual sources of loss
###

g_distances_scatter_heritab = ggplot(dataframe_0.00003, aes( x=(1-exp(-25*heritab))   , y=distance)) +
  geom_point() +
  xlab("Probability of heritable tape loss")+
  ylab("Normalised distance to true tree")+
  theme_classic()+ 
  facet_wrap(facets=~metric) + 
  theme(strip.text.x = element_text(size = 11),plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = c(0.5,0.5))


g_distances_scatter_dropout = ggplot(dataframe_0.00003, aes( x=(exp(-25*heritab)*droping)    , y=distance)) +
  geom_point() +
  xlab("Probability of tape dropout")+
  ylab("Normalised distance to true tree")+
  theme_classic()+ 
  facet_wrap(facets=~metric) + 
  theme(strip.text.x = element_text(size = 11),plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = c(0.5,0.5))

both_plots_2 <- plot_grid(g_distances_scatter_heritab,g_distances_scatter_dropout,ncol=1)

ggsave(plot = both_plots_2, "plots/supp_fig_24.pdf",
       width = 14.28, height = 14.28, units = "cm")



clock_rate_inference <- read.csv("supp_fig_lossy_data/clock_rate_inference_no_missing_model.csv")
death_rate_inference <- read.csv("supp_fig_lossy_data/death_rate_inference_no_missing_model.csv")
clock_rate_inference <- read.csv("supp_fig_lossy_data/clock_rate_inference_no_missing_model.csv")
insert_rate_inference <- read.csv("supp_fig_lossy_data/insert_rate_inference_no_missing_model.csv")
growth_rate_inference <- read.csv("supp_fig_lossy_data/growth_rate_inference_no_missing_model.csv")

selected_clock_rate_inference <- clock_rate_inference[which(clock_rate_inference$median != 0.0),]
selected_death_rate_inference <- death_rate_inference[which(death_rate_inference$median != 0.0),]
selected_clock_rate_inference <- clock_rate_inference[which(clock_rate_inference$median != 0.0),]
selected_growth_rate_inference <- growth_rate_inference[which(growth_rate_inference$median != 0.0),]
selected_insert_rate_inference <- insert_rate_inference[which(insert_rate_inference$median != 0.0),]

bias_growth_without_modelling <- (selected_growth_rate_inference$median-0.6)
bias_insert_without_modelling <- (selected_insert_rate_inference$median-selected_insert_rate_inference$true_value)
bias_birth_without_modelling <- (selected_clock_rate_inference$median-0.8)
bias_death_without_modelling <- (selected_death_rate_inference$median-0.2)
bias_clock_without_modelling <- (selected_clock_rate_inference$median-selected_clock_rate_inference$true_value)

hpd_clock_without_modelling <- (abs(selected_clock_rate_inference$hpd_upper-selected_clock_rate_inference$hpd_lower)/selected_clock_rate_inference$true_value)
hpd_birth_without_modelling <- (abs(selected_clock_rate_inference$hpd_upper-selected_clock_rate_inference$hpd_lower)/selected_clock_rate_inference$true_value)
hpd_death_without_modelling <- (abs(selected_death_rate_inference$hpd_upper-selected_death_rate_inference$hpd_lower)/selected_death_rate_inference$true_value)
hpd_growth_without_modelling <- (abs(selected_growth_rate_inference$hpd_upper-selected_growth_rate_inference$hpd_lower)/selected_growth_rate_inference$true_value)
hpd_insert_without_modelling <- (abs(selected_insert_rate_inference$hpd_upper-selected_insert_rate_inference$hpd_lower)/selected_insert_rate_inference$true_value)


all_stats_growth_editing_rate_without_loss <- data.frame(bias=c(bias_clock_without_modelling,bias_insert_without_modelling ,bias_growth_without_modelling),hpd=c(hpd_clock_without_modelling,hpd_insert_without_modelling,hpd_growth_without_modelling),parameter=c(rep("Editing rate",length(bias_clock_without_modelling)),rep("Insert probabilities",length(bias_insert_without_modelling)),rep("Growth rate",length(bias_growth_without_modelling)))) 





###
# Bias/uncertainty
###

g_bias = ggplot(all_stats_growth_editing_rate_without_loss, aes( x=factor(parameter), y=bias)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1,  alpha=0.2, outlier.shape = NA,coef = 0) +
  xlab("Parameter")+
  ylab("Bias of posterior median")+
  theme_classic()+ 
  ylim(-0.6,0.6) +
  ggtitle("Bias in estimates")+ 
  theme(strip.text.x = element_text(size = 11),plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = c(0.5,0.5))

g_hpd = ggplot(all_stats_growth_editing_rate_without_loss, aes( x=factor(parameter), y=hpd)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1,  alpha=0.2, outlier.shape = NA,coef = 0) +
  xlab("Parameter")+
  ylab("Relative HPD width")+
  theme_classic()+ 
  ylim(0.1,2.0) +
  ggtitle("Uncertainty in estimates") +
  theme(strip.text.x = element_text(size = 11),plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = c(0.5,0.5))

#combine in grid
both_plots_3 <- plot_grid(g_bias,g_hpd,nrow=1)

ggsave(plot = both_plots_3 + theme(text=element_text(size = text_size)), "plots/supp_fig_25.pdf",
       width = 14.28, height = 7.14, units = "cm")


g_coverage_all = ggplot(filtering_stats_0.00003, aes( x=mean_coverage)) +
  geom_histogram(bins = 40,col="#D1CD92",fill="grey") +
  xlab("Mean total tapes recovered per cell")+
  ylab("Number of datasets")+
  theme_classic()+ 
  theme(strip.text.x = element_text(size = 11),plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = c(0.5,0.9))

ggsave(plot = g_coverage_all + theme(text=element_text(size = text_size)), "plots/supp_fig_26.pdf",
       width = 14.28, height = 7.14, units = "cm")


