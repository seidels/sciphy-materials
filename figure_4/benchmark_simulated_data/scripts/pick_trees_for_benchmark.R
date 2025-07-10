# Choose 9 trees and alignments for benchmark
# Choose 3 trees for each tree size category, ranging from 
#.    small trees with <= 150 leaves
#.    medium trees with more than 150 and <= 450 leaves
#.    large trees with more than 450 leaves

# import libs
library(ape)
library(treeio)

# init path
path_to_trees = "benchmark/validations_2024/data/"

# init variables
set.seed(6)
tree_list = vector("list", 100)
tree_dat = data.frame(tree_number=1:100, nNodes=0)

# src function to read newick tree file that is missing a closing semicolon
source("./benchmark/scripts/functions.R")

# Collect all trees from validation study 
for (tree_number in 1:100){
  tree_file = paste0(path_to_trees, "simulate_alignment_and_tree.", tree_number, ".newick")
  tree = read_tree_with_fix(tree_file)
  tree_list[[tree_number]] = tree 
  tree_dat[tree_number, "nNodes"] = tree$Nnode
  print(tree_file)
}

# Choose 3 trees for each tree size category
## Note that nNodes is the number of internal nodes, and the number of leaves is
## thus nNodes + 1.
# The print statement shows the tree numbers that are chosen for each category.
trees_small = tree_dat[which(tree_dat$nNodes <= 149), ]
trees_small_selection = sample(trees_small$tree_number, size = 3)
trees_small_selection
print(tree_dat[which(tree_dat$tree_number %in% trees_small_selection), ]) # 13, 61, 72

trees_medium = tree_dat[which(tree_dat$nNodes > 149 & tree_dat$nNodes <= 449), ]
trees_medium_selection = sample(trees_medium$tree_number, size = 3)
trees_medium_selection
print(tree_dat[which(tree_dat$tree_number %in% trees_medium_selection), ]) # 15, 17, 60

trees_large = tree_dat[which(tree_dat$nNodes > 449), ]
trees_large_selection = sample(trees_large$tree_number, size = 3)
trees_large_selection
print(tree_dat[which(tree_dat$tree_number %in% trees_large_selection), ]) # 2, 22, 98
