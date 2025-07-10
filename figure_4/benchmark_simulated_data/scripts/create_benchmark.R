
library(TreeDist)
library(ape)
library(stringr)
source("benchmark/scripts/functions.R")
library(phangorn)

path_to_upgma_output = "benchmark/upgma-output/"
path_to_ordered_upgma_output = "benchmark/ordered-upgma-output/"

path_to_random_tree_output = "benchmark/random_tree_output/"

path_to_tidetree_output = "benchmark/tidetree-output/"
path_to_sciphy_output = "benchmark/validations_2024/log/"

path_to_data = "benchmark/alignments/"
path_to_sciphy_mcc = "benchmark/sciphy_output/MCC/"
path_to_sciphy_ccd = "benchmark/sciphy_output/CCD/"
tree_categories = list.dirs(path_to_upgma_output, recursive = F, full.names = F)
tree_categories

distance_results = data.frame(tree_categories = rep("", 9),
                              tree_number = 0, 
                              RF_distance_upgma = -1, 
                              RF_distance_tidetree_mcc = -1,
                              RF_distance_tidetree_ccd = -1,
                              RF_distance_sciphy_mcc = -1,
                              RF_distance_sciphy_ccd = -1,
                              RF_distance_upgma_ordered = -1,
                              RF_distance_random = -1,
                              wRF_distance_upgma = -1, 
                              wRF_distance_tidetree_mcc = -1,
                              wRF_distance_tidetree_ccd = -1,
                              wRF_distance_sciphy_mcc = -1,
                              wRF_distance_sciphy_ccd = -1,
                              wRF_distance_upgma_ordered = -1,
                              wRF_distance_random = -1, 
                              PI_distance_upgma = -1, 
                              PI_distance_tidetree_mcc = -1,
                              PI_distance_tidetree_ccd = -1,
                              PI_distance_sciphy_mcc = -1,
                              PI_distance_sciphy_ccd = -1,
                              PI_distance_upgma_ordered = -1,
                              PI_distance_random = -1)

ctr = 1
for (tree_category in tree_categories){
  #tree_category = tree_categories[3]
  print(tree_category)
  upgma_tree_files = list.files(paste0(path_to_upgma_output, tree_category))
  upgma_tree_files
  
  for (tree_file in upgma_tree_files){
    #tree_file = upgma_tree_files[1]
    print(tree_file)
    
    # get upgma tree
    upgma_tree_file = paste0(path_to_upgma_output, tree_category, "/", tree_file)
    print(upgma_tree_file)
    upgma_tree = read.tree(file = upgma_tree_file)
    
    tree_number = str_extract(tree_file, "(?<=upgma-)\\d+(?=\\.9\\.alignment\\.newick)")
    tree_number
    print(tree_number)
    
    # get ordered upgma tree
    ordered_upgma_tree_file = paste0(path_to_ordered_upgma_output, tree_category, "/upgma-ordered-", tree_number, ".newick")
    ordered_upgma_tree = read.tree(file = ordered_upgma_tree_file)
    
    # get random tree
    random_tree_file = paste0(path_to_random_tree_output, tree_category, "/bd-tree-", tree_number, ".newick")
    random_tree = read.tree(file = random_tree_file)
    
    # get tidetree tree
    tidetree_mcc = ape::read.nexus(file = paste0(path_to_tidetree_output, 
                                         "/infer_given_fixed_tree_13_inserts.tree.", tree_number, ".mcc.tree"))
    tidetree_ccd = ape::read.nexus(file = paste0(path_to_tidetree_output, 
                                                 "/infer_given_fixed_tree_13_inserts.tree.", tree_number, ".ccd0.tree"))
    
    # get sciphy tree
    sciphy_mcc = ape::read.nexus(file = paste0(path_to_sciphy_mcc, 
                                                  "MCC_tree.", tree_number, ".txt"))
    sciphy_ccd = ape::read.nexus(file = paste0(path_to_sciphy_ccd, 
                                               "CCD_tree.", tree_number, ".txt"))
    
    # get true tree
    true_tree_file = paste0(path_to_data, tree_category, "/simulate_alignment_and_tree.", tree_number, ".newick")   
    true_tree = read_tree_with_fix(true_tree_file)
    
    #plot_scaled_tree(tree, 25)
    
    ## Calculate distances
    #### UPGMA
    dist_upgma = TreeDist::RobinsonFoulds(tree1 = upgma_tree, tree2 = true_tree, normalize = T)
    wRF_upgma = wRF.dist(tree1 = upgma_tree, tree2 = true_tree, normalize = T, rooted = T)
    pi_upgma = TreeDist::PhylogeneticInfoDistance(tree1 = upgma_tree, tree2 = true_tree, normalize = T)
    distance_results[ctr, "tree_categories"] = tree_category
    distance_results[ctr, "tree_number"] = tree_number
    distance_results[ctr, "RF_distance_upgma"] = dist_upgma
    distance_results[ctr, "PI_distance_upgma"] = pi_upgma
    distance_results[ctr, "wRF_distance_upgma"] = wRF_upgma
    
    #### ordered UPGMA
    dist_upgma_ordered = TreeDist::RobinsonFoulds(tree1 = ordered_upgma_tree, tree2 = true_tree, normalize = T)
    pi_upgma_ordered = TreeDist::PhylogeneticInfoDistance(tree1 = ordered_upgma_tree, tree2 = true_tree, normalize = T)
    wRF_upgma_ordered = wRF.dist(tree1 = ordered_upgma_tree, tree2 = true_tree, normalize = T, rooted = T)
    distance_results[ctr, "RF_distance_upgma_ordered"] = dist_upgma_ordered
    distance_results[ctr, "wRF_distance_upgma_ordered"] = wRF_upgma_ordered
    distance_results[ctr, "PI_distance_upgma_ordered"] = pi_upgma_ordered
    
    #### random tree
    dist_random_tree = TreeDist::RobinsonFoulds(tree1 = random_tree, tree2 = true_tree, normalize = T)
    pi_random = TreeDist::PhylogeneticInfoDistance(tree1 = random_tree, tree2 = true_tree, normalize = T)
    wRF_random = wRF.dist(tree1 = random_tree, tree2 = true_tree, normalize = T, rooted = T, check.labels = F)
    
    distance_results[ctr, "RF_distance_random"] = dist_random_tree
    distance_results[ctr, "wRF_distance_random"] = wRF_random
    distance_results[ctr, "PI_distance_random"] = pi_random
    
    #### TiDeTree
    dist_tidetree_mcc = TreeDist::RobinsonFoulds(tree1 = tidetree_mcc, tree2 = true_tree, normalize = T)
    pi_tidetree_mcc = TreeDist::PhylogeneticInfoDistance(tree1 = tidetree_mcc, tree2 = true_tree, normalize = T)
    wRF_tidetree_mcc = wRF.dist(tree1 = tidetree_mcc, tree2 = true_tree, normalize = T, rooted = T)
    
    distance_results[ctr, "RF_distance_tidetree_mcc"] = dist_tidetree_mcc
    distance_results[ctr, "wRF_distance_tidetree_mcc"] = wRF_tidetree_mcc
    distance_results[ctr, "PI_distance_tidetree_mcc"] = pi_tidetree_mcc
    
    dist_tidetree_ccd = TreeDist::RobinsonFoulds(tree1 = tidetree_ccd, tree2 = true_tree, normalize = T)
    pi_tidetree_ccd = TreeDist::PhylogeneticInfoDistance(tree1 = tidetree_ccd, tree2 = true_tree, normalize = T)
    wRF_tidetree_ccd = wRF.dist(tree1 = tidetree_ccd, tree2 = true_tree, normalize = T, rooted = T)
    
    distance_results[ctr, "RF_distance_tidetree_ccd"] = dist_tidetree_ccd
    distance_results[ctr, "wRF_distance_tidetree_ccd"] = wRF_tidetree_ccd
    distance_results[ctr, "PI_distance_tidetree_ccd"] = pi_tidetree_ccd
    
    
    #### SciPhy
    dist_sciphy_mcc = TreeDist::RobinsonFoulds(tree1 = sciphy_mcc, tree2 = true_tree, normalize = T)
    pi_sciphy_mcc = TreeDist::PhylogeneticInfoDistance(tree1 = sciphy_mcc, tree2 = true_tree, normalize = T)
    wRF_sciphy_mcc = wRF.dist(tree1 = sciphy_mcc, tree2 = true_tree, normalize = T, rooted = T)
    
    
    distance_results[ctr, "RF_distance_sciphy_mcc"] = dist_sciphy_mcc
    distance_results[ctr, "wRF_distance_sciphy_mcc"] = wRF_sciphy_mcc
    distance_results[ctr, "PI_distance_sciphy_mcc"] = pi_sciphy_mcc
    
    dist_sciphy_ccd = TreeDist::RobinsonFoulds(tree1 = sciphy_ccd, tree2 = true_tree, normalize = T)
    pi_sciphy_ccd = TreeDist::PhylogeneticInfoDistance(tree1 = sciphy_ccd, tree2 = true_tree, normalize = T)
    wRF_sciphy_ccd = wRF.dist(tree1 = sciphy_ccd, tree2 = true_tree, normalize = T, rooted = T)
    
    
    distance_results[ctr, "RF_distance_sciphy_ccd"] = dist_sciphy_ccd
    distance_results[ctr, "wRF_distance_sciphy_ccd"] = wRF_sciphy_ccd
    distance_results[ctr, "PI_distance_sciphy_ccd"] = pi_sciphy_ccd
    
    
    ctr = ctr + 1
    
    }
}

saveRDS(distance_results, file = "benchmark/distances_across_methods.RDS")
