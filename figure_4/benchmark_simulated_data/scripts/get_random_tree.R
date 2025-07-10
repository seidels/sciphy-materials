library(ape)
library(stringr)
source("benchmark/scripts/functions.R")

path_to_data = "benchmark/alignments/"
path_to_output = "benchmark/random_tree_output/"
subdirs = list.dirs(path_to_data, full.names = F, recursive = F)
subdirs = subdirs[subdirs != "all"]
subdirs # = large_trees, medium_trees, small_trees

set.seed(1)
for (tree_category in subdirs){
  path_to_data_category = paste0(path_to_data, tree_category, "/")
  tree_files = list.files(path_to_data_category, pattern = "*.newick")
  tree_numbers = str_extract(tree_files, "(?<=tree\\.)\\d+(?=\\.newick)")
  
  tree_numbers
  for (tree_number in tree_numbers){
    
    # get true tree
    true_tree_file = paste0(path_to_data, tree_category, "/simulate_alignment_and_tree.", tree_number, ".newick")   
    true_tree = read_tree_with_fix(true_tree_file)
    
    n_tips = length(true_tree$tip.label)
    
    print("Getting random tree")
    random_tree <-  TreeSim::sim.bd.taxa.age(n_tips,1,0.8,0.2,0.00003,25)
    #adjusting tip labels such that the distance can be calculated
    random_tree[[1]]$tip.label <- str_remove(random_tree[[1]]$tip.label,"t")
    
    output_file = paste0(path_to_output,
                         tree_category,
                         "/bd-tree-",
                         tree_number,
                         ".newick"
    )
    write.tree(random_tree, file = output_file)
    
  }
}
