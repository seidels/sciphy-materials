library(ape)
library(stringr)
source("benchmark/scripts/functions.R")

path_to_data = "benchmark/alignments/"
path_to_output = "benchmark/upgma-output/"
subdirs = list.dirs(path_to_data, full.names = F, recursive = F)
subdirs = subdirs[subdirs != "all"]
subdirs # = large_trees, medium_trees, small_trees

for (tree_category in subdirs){
  print(paste0("Getting alignments from the ", tree_category))
  path_to_data_category = paste0(path_to_data, tree_category, "/")
  tree_files = list.files(path_to_data_category, pattern = "*.newick")
  tree_numbers = str_extract(tree_files, "(?<=tree\\.)\\d+(?=\\.newick)")
  
  tree_numbers
  for (tree_number in tree_numbers){
    
    alignment_files = list.files(path_to_data_category,  
                                 pattern = paste0("simulate_alignment_and_tree\\.seed=", tree_number, "\\.[0-9]+.alignment.nexus"))
  
    alignment_files  
    alignment_list <- list()
    for (alignment_file in alignment_files) {
      path_to_file <- file.path(path_to_data_category, alignment_file)
      print(paste("Reading in alignment file ", path_to_file))
      
      alignment <- read_integer_alignment(path_to_file)
      alignment_list[[alignment_file]] <- alignment
    }
    
    # Check that all alignments have identical row names
    if (!all(sapply(alignment_list, function(x) identical(rownames(x), rownames(alignment_list[[1]]))))) {
      stop("Row names (taxa) differ between alignments.")
    }
    
    # Concatenate all alignments by columns
    combined_alignment <- do.call(cbind, alignment_list)
    colnames(combined_alignment) = paste0("site_", 1:50)
    
    print("Computing distance matrix")
    dist_mat = get_distance_matrix(combined_alignment)
    dist_obj = as.dist(dist_mat)
    
    print("Getting upgma tree")
    upgma_tree = as.phylo(hclust(d = dist_obj, method = "average"))
    upgma_tree$edge.length[upgma_tree$edge.length == 0] <- 0.0001
    scaled_tree = scale_tree(tree = upgma_tree, root_heigh = 25)
    plot_scaled_tree(scaled_tree, 25)
    output_file = paste0(path_to_output,
                         tree_category,
                         "/upgma-",
                         strsplit(strsplit(alignment_file, "=")[[1]][2], ".nexus")[[1]][1],
                         ".newick"
    )
    write.tree(upgma_tree, file = output_file)
    
  }
}

#     
#   }
#   for (alignment_file in alignment_files){
#     path_to_file = paste0(path_to_data_category, alignment_file)
#     
#     print(paste0("Reading in alignment file ", path_to_file))
#     alignment = read_integer_alignment(path_to_file)
#     
#     print("Computing distance matrix")
#     dist_mat = get_distance_matrix(alignment)
#     dist_obj = as.dist(dist_mat)
#     
#     print("Getting upgma tree")
#     upgma_tree = as.phylo(hclust(d = dist_obj, method = "average"))
#     scaled_tree = scale_tree(tree = upgma_tree, root_heigh = 25)
#     plot_scaled_tree(scaled_tree, 25)
#     output_file = paste0(path_to_output,
#                          tree_category,
#                          "/upgma-",
#                          strsplit(strsplit(alignment_file, "=")[[1]][2], ".nexus")[[1]][1],
#                          ".newick"
#                          )
#     write.tree(scaled_tree, file = output_file)
#   }
# }
# 
# 
