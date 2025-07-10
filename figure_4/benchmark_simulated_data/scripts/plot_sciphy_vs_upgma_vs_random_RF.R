## ---------------------------
##
## Script name: Validation dataset features + benchmarking
##
## Purpose of script: explore features of the dataset simulated for the 100 alignments validation + benchmarking
##
## Author: Antoine Zwaans
##
## Date Created: 2024-03-05
##
## Copyright (c) Antoine Zwaans, 2024
## Email: antoine.zwaans@bsse.ethz.ch
##
## set working directory for the validation folder

setwd("/Users/azwaans/Documents/Projects/TYPEWRITER/validations_2025")   


## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
library(ape)
library(treebalance)
library(TreeSim)
library(rstatix)
library(ggpubr)
library(phylobase)

#load the upgma building function from Choi et al.,2022
source("~/typewriter_analysis/src/cell_culture/upgma/scripts.R")
source("~/typewriter_analysis/src/useful_scripts_across_categories.R")

#get num tips per simulated tree
num_tips <- c()
for(i in 1:100) {
  file_name <- paste0("../data_validations/simulate_alignment_and_tree.",i,".newick")
  tree <- readLines(file_name)
  tree <- str_remove_all(tree,"\\[.........\\]")
  tree <- str_remove_all(tree,"\\[..........\\]")
  tree <- str_remove_all(tree,"\\[...........\\]")
  tree <- paste0(tree,";")
  p_tree <- ape::read.tree(text = tree)
  num_tips <- c(num_tips,length(p_tree$tip.label))
}

#plot num tips per seed
num_tips <- data.frame(num=num_tips)
plot <- ggplot(data=num_tips,aes(x=num)) + geom_histogram(color="turquoise") + geom_vline(xintercept=median(num_tips$num)) + xlab("Number of tips") + ylab("Number of simulated trees") +theme_bw() + theme(text = element_text(size = 22)) 
ggsave("validation_2024_dataset_features.png", plot , width = 40, height =40 , units = "cm", dpi = 800)

## function to read-in the simulated alignments
parseSimulatedTapeAlignment <- function(seed) {
  full_data <- data.frame()
  for(tapeNR in 1:10) {
    file_name <- paste0("../data_validations/simulate_alignment_and_tree.seed=",seed,".",tapeNR,".alignment.nexus")
    tree <- readLines(file_name)
    tree <- tree[13:length(tree)-1]
    tree <- str_remove_all(tree,"\t\t")
    tree <- str_remove_all(tree,";")
    
    taxa <- c()
    sequences <- c()
    for(i in 1:length(tree)) {
      taxa <- c(taxa,unlist(str_split(tree[i]," "))[1]) 
      sequences <- c(sequences,unlist(str_split(tree[i]," "))[2]) 
    }
    
    sequences <- str_replace_all(sequences,"0","None")
    #putting back 10s that were removed. 
    sequences <- str_replace_all(sequences,"1None","10")
    sequences <- str_split(sequences,",",simplify = T)
    
    colnames(sequences) <- c("Site1","Site2","Site3","Site4","Site5")
    new_sequences <- data.frame(cbind(Cell=taxa,TargetBC=tapeNR,sequences))
    full_data <- rbind(full_data,new_sequences)
    
  }
  return(full_data)
}



### extract and plot distances to true tree: 
distance_UPGMA_truth_rf <- c()
distance_SCIPHY_truth_rf <- c()
distance_truth_random_rf <- c()

#distance_UPGMA_truth_PI <- c()
#distance_SCIPHY_truth_PI <- c()
#distance_truth_random_PI <- c()

for(SEED in 1:100) {
  
  simulated_alignment <- parseSimulatedTapeAlignment(seed=SEED)
  # nSites is 13 targetBCs * 5 sites = 65 sites; minus the contracted sites,
  # i.e. 1x 3 contracted (2xTape) and 3x 1 contracted (4xTape)
  #in our simulation setup, we have 50 sites
  UPGMA_tree = build_upgma_tree(edit_table = simulated_alignment, nSites = 50)
  
  #read-in the SciPhy MCC 
  MCC_tree <- ape::read.nexus(paste0("MCC/MCC_tree.",SEED,".txt"))
  
  #parse the true tree
  file_name_true_tree <- paste0("../data_validations/simulate_alignment_and_tree.",SEED,".newick")
  true_tree <- readLines(file_name_true_tree)
  true_tree <- str_remove_all(true_tree,"\\[.........\\]")
  true_tree <- str_remove_all(true_tree,"\\[..........\\]")
  true_tree <- str_remove_all(true_tree,"\\[...........\\]")
  true_tree <- paste0(true_tree,";")
  true_tree <- ape::read.tree(text = true_tree)
  
  #random tree simulated with the same b-d-sampling parameters as for the validation + same number of tips
  random_tree <-  TreeSim::sim.bd.taxa.age(length(true_tree$tip.label),1,0.8,0.2,0.00003,25)
  #adjusting tip labels such that the distance can be calculated
  random_tree[[1]]$tip.label <- str_remove(random_tree[[1]]$tip.label,"t")
  
  #calculate pairwise RF distances for each 
  distance_UPGMA_truth_rf <-  c(distance_UPGMA_truth_rf,TreeDist::RobinsonFoulds(true_tree, UPGMA_tree,normalize = TRUE))
  distance_SCIPHY_truth_rf <- c(distance_SCIPHY_truth_rf,TreeDist::RobinsonFoulds(true_tree, MCC_tree,normalize = TRUE))
  distance_truth_random_rf <- c(distance_truth_random_rf,TreeDist::RobinsonFoulds(random_tree, true_tree,normalize = TRUE))
  
  #calculate pairwise PI distances for each 
  #distance_UPGMA_truth <-  c(distance_UPGMA_truth,TreeDist::PhylogeneticInfoDistance(true_tree, seed1_UPGMA,normalize = TRUE))
  #distance_SCIPHY_truth <- c(distance_SCIPHY_truth,TreeDist::PhylogeneticInfoDistance(true_tree, seed1_MCC,normalize = TRUE))
  #distance_truth_random <- c(distance_truth_random,TreeDist::PhylogeneticInfoDistance(random_tree, true_tree,normalize = TRUE))
}

#make dataframe for boxplot + t-test
rf_distances_to_truth <- data.frame(seed=1:100,sciphy= distance_SCIPHY_truth_rf,upgma=distance_UPGMA_truth_rf,random=distance_truth_random_rf)
rf_distances_to_truth <- pivot_longer(rf_distances_to_truth, !seed,names_to = "reference", values_to = "distance")

#pariwise t-test
stat.test <- pairwise_t_test(distance ~ reference, data = rf_distances_to_truth, paired = TRUE, 
                             p.adjust.method = "bonferroni"
)
stat.test <- stat.test %>% add_xy_position(x = "reference")

bxp <- ggboxplot(
  rf_distances_to_truth, x = "reference", y = "distance",
  color = "reference", palette = "jco"
)+ ylab("RF distance to true tree") + xlab("Tree reconstruction method")# Add pairwise comparisons p-value
paired_test <- bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08
)
paired_test
ggsave("test_mean_distances_RF_paired_MCC.pdf",paired_test, width = 12, height = 12, units = "cm")