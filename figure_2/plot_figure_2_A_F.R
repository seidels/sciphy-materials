## ---------------------------
##
## Script name: get plot for all inferred parameters in the full validation
##
## Purpose of script: plot well-calibrated simulations
##
## Author: Antoine Zwaans, Sophie Seidel
##
## Date Created: 2024-07-05
##
## Copyright (c) Antoine Zwaans, Sophie Seidel, 2024
## Email: azwaans@ethz.ch, ophie.seidel@posteo.de
##
##
## ---------------------------

#text size for plotting
text_size = 10

#load packages
require(data.table)
library(ape)
library(treebalance)
library(TreeSim)
library(rstatix)
library(ggpubr)
library(phylobase)
library(tracerer)
library(HDInterval)
library(ggplot2)
library(beastio)
library(stringr)
library(tidyverse)
library(phytools)

#read-in validation results for each SciPhy parameter and tree statistic
B1_index_inference <- read.csv(file = "inference_logs/summary/B1_index_inference.csv")
insert_rate_inference <- read.csv(file = "inference_logs/summary/insert_rate_inference.csv")
clock_rate_inference <- read.csv(file = "inference_logs/summary/clock_rate_inference.csv")
tree_height_inference <- read.csv(file = "inference_logs/summary/tree_height_inference.csv")
tree_length_inference <- read.csv(file = "inference_logs/summary/tree_length_inference.csv")

# all chains converged, assessed with get_summary_tables_for_all_parameters.R
nr_converged_chains <- 100

#coverages per insert rate
coverages_per_insert <- c()
for(i in 1:13) {
  coverage <- sum(insert_rate_inference[which(insert_rate_inference$insertRate == i),"recovered"])/nr_converged_chains
  coverages_per_insert <- c(coverages_per_insert,coverage)
}
coverages_per_insert
#correlations per insert rate
for(i in 1:13) {
  print(cor.test(insert_rate_inference[which(insert_rate_inference$insertRate == i),"median"], insert_rate_inference[which(insert_rate_inference$insertRate == i),"true_value"], method = "pearson"))
}
#coverage for the clock rate
sum(clock_rate_inference$recovered)/nr_converged_chains
#correlation for the clock rate
correlation_clock <- cor.test(clock_rate_inference$median, clock_rate_inference$true_value, method = "pearson")

#coverage for the tree height
sum(tree_height_inference$recovered)/nr_converged_chains
#correlation for the tree height
correlation_tree_height <- cor.test(tree_height_inference$median, tree_height_inference$true_value, method = "pearson")

#coverage for the tree length
sum(tree_length_inference$recovered)/nr_converged_chains
#correlation for the tree length
correlation_tree_length <- cor.test(tree_length_inference$median,tree_length_inference$true_value, method = "pearson")

#coverage for the B1 index
sum(B1_index_inference$recovered)/nr_converged_chains
#correlation for the B1 index
correlation_B1 <- cor.test(B1_index_inference$median,B1_index_inference$true_value, method = "pearson")

## true values on x axis and estimated values on y axis
cols <- c("TRUE" = "black", "FALSE" = "darkgrey")

insert_probs_plot = ggplot(insert_rate_inference, aes(x=true_value, y=median,color=as.logical(recovered))) +
  geom_point(size=0.3)+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper,color=as.logical(recovered)), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median \n& posterior interval")+
  theme_classic()+
  xlim(0, 0.5) + ylim(0, 0.5)+
  ggtitle("Insertion probabilities")+ 
  scale_color_manual(values = cols) +
  theme(plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = "None")
#insert_probs_plot

clock_rate_plot = ggplot(clock_rate_inference, aes(x=true_value, y=median,color=as.logical(recovered))) +
  geom_point(size=0.3)+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper,color=as.logical(recovered)), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median \n& posterior interval")+
  theme_classic()+
  xlim(0, 0.5) + ylim(0, 0.5)+
  ggtitle("Editing rate")+
  scale_color_manual(values = cols) +
  theme(plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = "None")
#clock_rate_plot

tree_height_plot = ggplot(tree_height_inference, aes(x=true_value, y=median,color=as.logical(recovered))) +
  geom_point(size=0.3)+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper,color=as.logical(recovered)), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median \n& posterior interval")+
  theme_classic()+
  coord_cartesian(ylim=c(15, 30),xlim=c(15, 30)) +
  ggtitle("Tree height")+ 
  scale_color_manual(values = cols) +
  theme(plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = "None")
#tree_height_plot

tree_length_plot = ggplot(tree_length_inference, aes(x=true_value, y=median,color=as.logical(recovered))) +
  geom_point(size=0.3)+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper,color=as.logical(recovered)), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median\n& posterior interval")+
  theme_classic()+
  coord_cartesian(ylim=c(0, 10000),xlim=c(0, 10000)) +
  ggtitle("Tree length")+ 
  scale_color_manual(values = cols) +
  theme(plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = "None")

#tree_length_plot
B1_index_plot = ggplot(B1_index_inference, aes(x=true_value, y=median, color=as.logical(recovered))) +
  geom_point(size=0.3)+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper,color=as.logical(recovered)), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median \n& posterior interval")+
  theme_classic()+ 
  coord_cartesian(ylim=c(0,220 ),xlim=c(0, 220)) +
  ggtitle("Tree balance")+ 
  scale_color_manual(values = cols) +
  theme(plot.title = element_text(hjust=0.5),text=element_text(size = text_size),legend.position = "None")
#B1_index_plot

#get num tips per simulated tree
num_tips <- c()
for(i in 1:100) {
  file_name <- paste0("simulated_data/simulate_alignment_and_tree.",i,".newick")
  tree <- readLines(file_name)
  tree <- str_remove_all(tree,"\\[.........\\]")
  tree <- str_remove_all(tree,"\\[..........\\]")
  tree <- str_remove_all(tree,"\\[...........\\]")
  tree <- paste0(tree,";")
  p_tree <- ape::read.tree(text = tree)
  num_tips <- c(num_tips,length(p_tree$tip.label))
}


### extract and plot distances from CCD summary tree for each inference run to true tree in PI
distance_SCIPHY_truth_PI <- c()
distance_truth_random_PI <- c()

for(SEED in 1:100) {
  
  #read-in the SciPhy CCD
  CCD_tree <- ape::read.nexus(paste0("inference_logs/CCD/CCD_tree.",SEED,".txt"))
  
  #parse the true simulated tree
  file_name_true_tree <- paste0("simulated_data/simulate_alignment_and_tree.",SEED,".newick")
  true_tree <- readLines(file_name_true_tree)
  true_tree <- str_remove_all(true_tree,"\\[.........\\]")
  true_tree <- str_remove_all(true_tree,"\\[..........\\]")
  true_tree <- str_remove_all(true_tree,"\\[...........\\]")
  true_tree <- paste0(true_tree,";")
  true_tree <- ape::read.tree(text = true_tree)
  
  #Simulate and random tree with the same b-d-sampling parameters as for the validation + same number of tips
  random_tree <-  TreeSim::sim.bd.taxa.age(length(true_tree$tip.label),1,0.8,0.2,0.00003,25)
  
  #adjusting tip labels such that the distance can be calculated
  random_tree[[1]]$tip.label <- str_remove(random_tree[[1]]$tip.label,"t")
  
  #calculate pairwise PI distances to true tree topology for both random and inferred tree
  distance_SCIPHY_truth_PI <- c(distance_SCIPHY_truth_PI,TreeDist::PhylogeneticInfoDistance(true_tree, CCD_tree,normalize = TRUE))
  distance_truth_random_PI <- c(distance_truth_random_PI,TreeDist::PhylogeneticInfoDistance(random_tree, true_tree,normalize = TRUE))
}

#make dataframe for boxplot + t-test
pi_distances_to_truth <- data.frame(seed=1:100,SciPhy_CCD= distance_SCIPHY_truth_PI,Random_BD_tree=distance_truth_random_PI)
pi_distances_to_truth <- pivot_longer(pi_distances_to_truth, !seed,names_to = "reference", values_to = "distance")

#pariwise t-test
stat.test <- pairwise_t_test(distance ~ reference, data = pi_distances_to_truth, paired = TRUE, 
                             p.adjust.method = "bonferroni"
)
stat.test <- stat.test %>% add_xy_position(x = "reference")

#specify colors
cols <- c("SciPhy_CCD" = "#2e7d32cc","Random_BD_tree" = "#005bf2cc")

bxp_CCD_PI <- ggboxplot(
  pi_distances_to_truth, x = "reference", y = "distance",
  fill = "reference", palette = "jco", size = 0.2, outlier.size = 0.3
)+ ylab("Normalised PI \ndistance to true tree") + xlab("Tree reconstruction method")# Add pairwise comparisons p-value
paired_test_CCD_PI <- bxp_CCD_PI + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.2,y.position = 1.1
) +theme_classic()+ theme(legend.position = "None",text=element_text(size = text_size),plot.title = element_text(hjust=0.5)) + 
  scale_fill_manual(values = cols)  +  ggtitle("Tree topology") + 
  scale_y_continuous(breaks=c(0.0,0.2,0.6,0.8,1.0)) + ylim(c(0.0,1.5))
paired_test_CCD_PI
full_figure <- cowplot::plot_grid(insert_probs_plot,clock_rate_plot,tree_height_plot,tree_length_plot,B1_index_plot,paired_test_CCD_PI, ncol=2, labels = "AUTO") 

ggsave("figure_2_A_F.pdf",full_figure, width = 14.28, height = 14.28, units = "cm", dpi = 300)
