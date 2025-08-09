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

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
library(ape)
library(treebalance)
library(TreeSim)
library(rstatix)
library(ggpubr)
library(phylobase)
library(phangorn)

#load the upgma building function from Choi et al.,2022
source("~/typewriter_analysis/src/cell_culture/upgma/scripts.R")
source("~/typewriter_analysis/src/useful_scripts_across_categories.R")

#get num tips per simulated tree
num_tips <- c()
collision_probs <- c()
clock_rates <- c()
true_heights <- c()
for(i in 1:100) {
  file_name <- paste0("../figure_2/simulated_data/simulate_alignment_and_tree.",i,".newick")
  tree <- readLines(file_name)
  tree <- str_remove_all(tree,"\\[.........\\]")
  tree <- str_remove_all(tree,"\\[..........\\]")
  tree <- str_remove_all(tree,"\\[...........\\]")
  tree <- paste0(tree,";")
  p_tree <- ape::read.tree(text = tree)
  num_tips <- c(num_tips,length(p_tree$tip.label))
  true_heights <- c(true_heights,max(nodeHeights(p_tree)))
  
  file_name <-  paste0("../figure_2/simulation_parameters/simParams_",i,".csv")
  parameters <- read.csv(file_name)
  collision_probs <- c(collision_probs,sum(parameters$x[1:13]^2))
  clock_rates <- c(clock_rates,parameters$x[14])
  
}

## function to read-in the simulated alignments
parseSimulatedTapeAlignment <- function(seed) {
  full_data <- data.frame()
  for(tapeNR in 1:10) {
    file_name <- paste0(".../figure_2/simulated_data/simulate_alignment_and_tree.seed=",seed,".",tapeNR,".alignment.sciphy")
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



build_unordered_UPGMA <- function(edit_table = simulated_alignment) {

  distance_matrix <- matrix(0,nrow=length(unique(simulated_alignment$Cell)),ncol=length(unique(simulated_alignment$Cell)))
  
  for(i in 1:length(unique(simulated_alignment$Cell))) { 
    print(i)
    for(j in i:length(unique(simulated_alignment$Cell))) {
      
      seq_cell_i <- c()
      seq_cell_j <- c()
      total_distance <- 0
      for(targ in unique(simulated_alignment$TargetBC)) { 
        cell_index_i <- as.character(i-1) 
        cell_index_j <- as.character(j-1)
       #concatenate: 
        seq_cell_i <- c(seq_cell_i,unlist(simulated_alignment[simulated_alignment$Cell == cell_index_i & simulated_alignment$TargetBC == targ  ,3:7]))
        seq_cell_j <- c(seq_cell_j,unlist(simulated_alignment[simulated_alignment$Cell == cell_index_j & simulated_alignment$TargetBC == targ  ,3:7]))
        
      }
      total_distance <- length(seq_cell_i[seq_cell_i != seq_cell_j])
      distance_matrix[i,j] <- total_distance
      distance_matrix[j,i] <- total_distance
    }

    
    
  }
  

  dist_obj = as.dist(distance_matrix)
  print("Getting upgma tree")
  unordered_UPGMA = as.phylo(hclust(d = dist_obj, method = "average"))
  unordered_UPGMA$tip.label <- as.character(as.numeric(unordered_UPGMA$tip.label) - 1)
  return(unordered_UPGMA)
  
  
}


##
# Scaling to 25 days
##


distance_UPGMA_truth_PI <- c()
distance_SCIPHY_truth_PI <- c()
distance_truth_random_PI <- c()
distance_truth_unordered_UPGMA_PI <- c()

distance_UPGMA_truth_wrf <- c()
distance_SCIPHY_truth_wrf <- c()
distance_truth_random_wrf <- c()
distance_truth_unordered_UPGMA_wrf <- c()

original_heigths_UPGMA_ordered <- c()
original_heigths_UPGMA_unordered <- c()


for(SEED in 1:100) {
  
  
  print(paste0("seed: ",SEED))
  
  simulated_alignment <- parseSimulatedTapeAlignment(seed=SEED)
  # nSites is 13 targetBCs * 5 sites = 65 sites; minus the contracted sites,
  # i.e. 1x 3 contracted (2xTape) and 3x 1 contracted (4xTape)
  #in our simulation setup, we have 50 sites
  ordered_UPGMA = build_upgma_tree(edit_table = simulated_alignment, nSites = 50)
  unordered_UPGMA = build_unordered_UPGMA(edit_table = simulated_alignment)
 
   #read-in the SciPhy CCD 
  CCD_tree <- ape::read.nexus(paste0("../figure_2/inference_logs/CCD/CCD_tree.",SEED,".txt"))

  #parse the true tree
  file_name_true_tree <- paste0("../figure_2/simulated_data/simulate_alignment_and_tree.",SEED,".newick")
  true_tree <- readLines(file_name_true_tree)
  true_tree <- str_remove_all(true_tree,"\\[.........\\]")
  true_tree <- str_remove_all(true_tree,"\\[..........\\]")
  true_tree <- str_remove_all(true_tree,"\\[...........\\]")
  true_tree <- paste0(true_tree,";")
  true_tree <- ape::read.tree(text = true_tree)
  
  
  original_heigths_UPGMA_ordered <- c(original_heigths_UPGMA_ordered,max(nodeHeights(ordered_UPGMA)))
  original_heigths_UPGMA_unordered <- c(original_heigths_UPGMA_unordered,max(nodeHeights(unordered_UPGMA)))
  
  #rescaling to 25
  unordered_UPGMA$edge.length <- unordered_UPGMA$edge.length * 25 / max(nodeHeights(unordered_UPGMA))
  ordered_UPGMA$edge.length <- ordered_UPGMA$edge.length *  25 / max(nodeHeights(ordered_UPGMA))
  
  
  
  #random tree simulated with the same b-d-sampling parameters as for the validation + same number of tips
  random_tree <-  TreeSim::sim.bd.taxa.age(length(true_tree$tip.label),1,0.8,0.2,0.00003,25)[[1]]
  #adjusting tip labels such that the distance can be calculated
  random_tree$tip.label <- str_remove(random_tree$tip.label,"t")
  random_tree$tip.label <- as.character(as.numeric(random_tree$tip.label) -1)
  
  newick_random <- unlist(write.tree(random_tree))
  random_formatted <- ape::read.tree(text = newick_random)
  
  #calculate pairwise RF distances for each
    distance_UPGMA_truth_PI <-  c(distance_UPGMA_truth_PI,TreeDist::PhylogeneticInfoDistance(true_tree, ordered_UPGMA,normalize = TRUE))
    distance_SCIPHY_truth_PI <- c(distance_SCIPHY_truth_PI,TreeDist::PhylogeneticInfoDistance(true_tree, CCD_tree,normalize = TRUE))
    distance_truth_random_PI <- c(distance_truth_random_PI,TreeDist::PhylogeneticInfoDistance( true_tree,random_formatted,normalize = TRUE))
    distance_truth_unordered_UPGMA_PI <- c(distance_truth_unordered_UPGMA_PI,TreeDist::PhylogeneticInfoDistance(true_tree,unordered_UPGMA,normalize = TRUE))
   
    

    distance_UPGMA_truth_wrf <-  c(distance_UPGMA_truth_wrf,wRF.dist(true_tree, ordered_UPGMA,normalize = TRUE))
    distance_SCIPHY_truth_wrf <- c(distance_SCIPHY_truth_wrf,wRF.dist(true_tree, CCD_tree,normalize = TRUE))
    distance_truth_random_wrf <- c(distance_truth_random_wrf,wRF.dist(true_tree, random_formatted,normalize = TRUE))
    distance_truth_unordered_UPGMA_wrf <- c(distance_truth_unordered_UPGMA_wrf,wRF.dist(true_tree,unordered_UPGMA,normalize = TRUE))



  
}

collision_categories <- rep("LOW",100)
collision_quantile <- quantile(collision_probs,probs=c(0.33333,0.666667))

collision_categories[which(collision_probs < collision_quantile[1])] <- "[0.09,0.11]"
collision_categories[which((collision_probs > collision_quantile[1]) & (collision_probs < collision_quantile[2]))] <- "[0.11,0.13]"
collision_categories[which( collision_probs > collision_quantile[2])] <-  "[0.13,0.21]"

clock_categories <- rep("LOW",100)
clock_quantile <- quantile(clock_rates,probs=c(0.33333,0.666667))

clock_categories[which(clock_rates <= clock_quantile[1])] <- "[0.04,0.10]"
clock_categories[which((clock_rates > clock_quantile[1]) & (clock_rates < clock_quantile[2]))] <- "[0.10,0.15]"
clock_categories[which( clock_rates >= clock_quantile[2])] <- "[0.15,0.38]"

tree_categories <- rep("LOW",100)
tree_quantile <- quantile(num_tips,probs=c(0.33333,0.666667))

tree_categories[which(num_tips <= tree_quantile[1])] <- "[0,50]"
tree_categories[which((num_tips > tree_quantile[1]) & (num_tips < tree_quantile[2]))] <- "[50,127]"
tree_categories[which( num_tips >= tree_quantile[2])] <- "[127,643]"


#make dataframe for boxplot + t-test
PI_distances_to_truth <- data.frame(seed=1:100,distance_SciPhy_CCD= distance_SCIPHY_truth_PI,distance_UPGMA_ordered=distance_UPGMA_truth_PI,collision_bins=collision_categories,clock_bins=clock_categories,tree_bins=tree_categories,distance_UPGMA_unordered=distance_truth_unordered_UPGMA_PI,distance_Random=distance_truth_random_PI)
PI_distances_to_truth <- pivot_longer(PI_distances_to_truth,
  cols = starts_with("distance"),
  names_to = "reference",
  names_prefix = "distance_",
  values_to = "distance",
  values_drop_na = TRUE
)

colors_fill <- c("SciPhy_CCD" = "#2e7d32cc","UPGMA_ordered" = "#ff7373ff","UPGMA_unordered" = "#e1af00cc","Random" = "#337bf3ff")

#pariwise t-test
stat.test <- pairwise_t_test(distance ~ reference, data = PI_distances_to_truth, paired = TRUE, 
                             p.adjust.method = "bonferroni"
)
stat.test <- stat.test %>% add_xy_position(x = "reference")

PI_distances_to_truth <- mutate(PI_distances_to_truth,collision_bins = fct_relevel(collision_bins, 
                          "[0.09,0.11]", "[0.11,0.13]", "[0.13,0.21]"))

PI_distances_to_truth <- mutate(PI_distances_to_truth,reference = fct_relevel(reference, 
                                                                                   "SciPhy_CCD", "UPGMA_ordered", "UPGMA_unordered","Random"))
bxp4 <- ggboxplot(
  PI_distances_to_truth, x = "collision_bins", y = "distance", fill="reference", palette = "jco"
)+ ylab("PI distance \nto true tree") + xlab("Collision probability bin") +
  scale_fill_manual(values=colors_fill)

# Add pairwise comparisons p-value
bxp4

PI_distances_to_truth <- mutate(PI_distances_to_truth,clock_bins = fct_relevel(clock_bins, 
                                                                               "[0.04,0.10]", "[0.10,0.15]", "[0.15,0.38]"))
bxp5 <- ggboxplot(
  PI_distances_to_truth, x = "clock_bins", y = "distance", fill="reference", palette = "jco"
)+ ylab("PI distance \nto true tree") + xlab("Editing rate bin")+
  scale_fill_manual(values=colors_fill)# Add pairwise comparisons p-value
bxp5
PI_distances_to_truth <- mutate(PI_distances_to_truth,tree_bins = fct_relevel(tree_bins, 
                                                                               "[0,50]", "[50,127]", "[127,643]"))
bxp6 <- ggboxplot(
  PI_distances_to_truth, x = "tree_bins", y = "distance", fill="reference", palette = "jco"
)+ ylab("PI distance \nto true tree") + xlab("Tree size bin") +
  scale_fill_manual(values=colors_fill)# Add pairwise comparisons p-value
bxp6



wrf_distances_to_truth <- data.frame(seed=1:100,distance_SciPhy_CCD= distance_SCIPHY_truth_wrf,distance_UPGMA_ordered=distance_UPGMA_truth_wrf,collision_bins=collision_categories,clock_bins=clock_categories,tree_bins=tree_categories,distance_UPGMA_unordered=distance_truth_unordered_UPGMA_wrf,distance_Random=distance_truth_random_wrf)
wrf_distances_to_truth <- pivot_longer(wrf_distances_to_truth,
                                      cols = starts_with("distance"),
                                      names_to = "reference",
                                      names_prefix = "distance_",
                                      values_to = "distance",
                                      values_drop_na = TRUE
)


wrf_distances_to_truth <- mutate(wrf_distances_to_truth,collision_bins = fct_relevel(collision_bins, 
                                                                                     "[0.09,0.11]", "[0.11,0.13]", "[0.13,0.21]"))
wrf_distances_to_truth <- mutate(wrf_distances_to_truth,reference = fct_relevel(reference, 
                                                                              "SciPhy_CCD", "UPGMA_ordered", "UPGMA_unordered","Random"))

bxp1 <- ggboxplot(
  wrf_distances_to_truth, x = "collision_bins", y = "distance", fill="reference", palette = "jco"
)+ ylab("wRF distance \nto true tree") + xlab("Collision probability bin") +
  scale_fill_manual(values=colors_fill)# Add pairwise comparisons p-value
bxp1

wrf_distances_to_truth <- mutate(wrf_distances_to_truth,clock_bins = fct_relevel(clock_bins, 
                                                                                 "[0.04,0.10]", "[0.10,0.15]", "[0.15,0.38]"))
bxp2 <- ggboxplot(
  wrf_distances_to_truth, x = "clock_bins", y = "distance", fill="reference", palette = "jco"
)+ ylab("wRF distance \nto true tree") + xlab("Editing rate bin") +
  scale_fill_manual(values=colors_fill)# Add pairwise comparisons p-value
bxp2
wrf_distances_to_truth <- mutate(wrf_distances_to_truth,tree_bins = fct_relevel(tree_bins, 
                                                                              "[0,50]", "[50,127]", "[127,643]"))
bxp3 <- ggboxplot(
  wrf_distances_to_truth, x = "tree_bins", y = "distance", fill="reference", palette = "jco"
)+ ylab("wRF distance \nto true tree") + xlab("Tree size bin") +
  scale_fill_manual(values=colors_fill) # Add pairwise comparisons p-value
bxp3

full_figure <- cowplot::plot_grid(bxp1 + theme(legend.title = element_blank(),legend.position = c(0.5,0.5),axis.text = element_text(size=9))+ guides(fill = guide_legend(ncol = 1)),bxp2+ theme(legend.position = "None",legend.title = element_blank(),axis.text = element_text(size=9)),bxp3+ theme(legend.position = "None",axis.text = element_text(size=9)),
                                  bxp4 + theme(legend.position = "None",legend.title = element_blank(),axis.text = element_text(size=9)),bxp5+ theme(legend.position = "None",legend.title = element_blank(),axis.text = element_text(size=9)),bxp6+ theme(legend.position = "None",axis.text = element_text(size=9)),
                                  nrow=2, labels = "AUTO") 

ggsave("supp_fig_accuracy_benchmark_UPGMA.pdf",full_figure, width = 21.42, height = 14.28, units = "cm", dpi = 300)





wrf_distances_sciphy_to_truth <- data.frame(seed=1:100,
                                            distance_SciPhy_CCD= distance_SCIPHY_truth_wrf,
                                            clocks=clock_rates,
                                            collisions=collision_probs,
                                            ntips=num_tips
                                            )

PI_distances_sciphy_to_truth <- data.frame(seed=1:100,
                                            distance_SciPhy_CCD= distance_SCIPHY_truth_PI,
                                            clocks=clock_rates,
                                            collisions=collision_probs,
                                            ntips=num_tips
)


##
# wRF space
##

corr_clock_wRF <- ggscatter(wrf_distances_sciphy_to_truth, x = "clocks", y = "distance_SciPhy_CCD",
                           #add = "reg.line",  # Add regressin line
                           #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                           #conf.int = TRUE # Add confidence interval
)  + ylab("SciPhy CCD \nwRF distance to truth") + xlab("Editing rates [d-1]") + theme_classic()


corr_coll_wRF <- ggscatter(wrf_distances_sciphy_to_truth, y = "distance_SciPhy_CCD", x = "collisions"
                          #add = "reg.line" , # Add regressin line
                          #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                          #conf.int = TRUE # Add confidence interval
)  + ylab("SciPhy CCD \neRF distance to truth") + xlab("Collision probability") + theme_classic()


corr_ntips_wRF <- ggscatter(wrf_distances_sciphy_to_truth, y = "distance_SciPhy_CCD", x = "ntips"
                           #add = "reg.line"  # Add regressin line
                           #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                           #conf.int = TRUE # Add confidence interval
)  + ylab("SciPhy CCD \nwRF distance to truth") + xlab("Tree size") + theme_classic()


##
# PI space
##


corr_clock_PI <- ggscatter(PI_distances_sciphy_to_truth, y = "distance_SciPhy_CCD", x = "clocks"
                       # add = "reg.line",  # Add regressin line
                        #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                        #conf.int = TRUE # Add confidence interval
)  + ylab("SciPhy CCD \nPI distance to truth") + xlab("Editing rates [d-1]") + theme_classic()


corr_coll_PI <- ggscatter(PI_distances_sciphy_to_truth, y = "distance_SciPhy_CCD", x = "collisions"
                      # add = "reg.line",  # Add regressin line
                      # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                      # conf.int = TRUE # Add confidence interval
)  + ylab("SciPhy CCD \nPI distance to truth") + xlab("Collision probability") + theme_classic()


corr_ntips_PI <- ggscatter(PI_distances_sciphy_to_truth, y = "distance_SciPhy_CCD", x = "ntips"
                        #add = "reg.line",  # Add regressin line
                       # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                        #conf.int = TRUE # Add confidence interval
)  + ylab("SciPhy CCD \nPI distance to truth") + xlab("Tree size") + theme_classic()


full_figure <- cowplot::plot_grid(corr_clock_wRF,corr_coll_wRF,corr_ntips_wRF,corr_clock_PI,corr_coll_PI,corr_ntips_PI,
                                  nrow=2, labels = "AUTO") 

ggsave("supp_fig_accuracy.pdf",full_figure, width = 21.42, height = 14.28, units = "cm", dpi = 300)


