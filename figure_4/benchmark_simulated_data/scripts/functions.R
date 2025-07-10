library(treeio)
library(tidyr)
library(dplyr)

read_tree_with_fix <- function(filepath) {
  # function to read newick tree file that is missing a closing semicolon
  # Read the tree as a single text string
  tree_text <- paste(readLines(filepath), collapse = "")
  
  # Add semicolon if missing
  if (!grepl(";$", tree_text)) {
    tree_text <- paste0(tree_text, ";")
  }
  
  # Read the tree from the corrected text
  read.newick(text = tree_text)
}


read_integer_alignment <- function(file) {
  lines <- readLines(file)
  
  # Find the start and end of the matrix
  matrix_start <- grep("^\\s*matrix", lines, ignore.case = TRUE)
  matrix_end <- grep("^\\s*end;", lines, ignore.case = TRUE)
  matrix_lines <- lines[(matrix_start + 1):(matrix_end[2] - 1)]
  
  # Remove empty lines
  matrix_lines <- matrix_lines[nzchar(matrix_lines)]
  
  parsed <- lapply(matrix_lines, function(line) {
    line <- trimws(line)
    line <- sub(";$", "", line)  # remove trailing semicolon
    
    if (!grepl("^\\S+\\s+\\S+", line)) return(NULL)  # skip malformed lines
    
    split_line <- strsplit(line, "\\s+")[[1]]
    id <- split_line[1]
    values <- as.integer(unlist(strsplit(split_line[2], ",")))
    list(id = id, values = values)
  })
  # Filter NULLs
  parsed <- Filter(Negate(is.null), parsed)
  
  taxon_ids <- sapply(parsed, `[[`, "id")
  value_matrix <- do.call(rbind, lapply(parsed, `[[`, "values"))
  
  alignment_df <- as.data.frame(value_matrix)
  rownames(alignment_df) <- make.unique(taxon_ids)
  
  return(alignment_df)
}

edit_distance <- function(x, y) {
  
  # Ensure same length
  stopifnot(length(x) == length(y))
  
  dist <- numeric(length(x))
  
  # Case 1: same values
  dist[x == y] <- 0
  
  # Case 2: one is 0, one is > 0
  dist[(x == 0 & y > 0) | (x > 0 & y == 0)] <- 1
  
  # Case 3: both > 0 and different
  dist[(x > 0 & y > 0) & (x != y)] <- 2
  
  sum(dist)
}

get_distance_matrix <- function(alignment){
  # Number of rows
  n <- nrow(alignment)
  
  # Initialize empty matrix for pairwise distances
  dist_matrix <- matrix(0, n, n)
  
  # Fill the matrix with Hamming distances
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist <- edit_distance(alignment[i, ], alignment[j, ])
      dist_matrix[i, j] <- dist
      dist_matrix[j, i] <- dist
    } 
  }
  
  # Set row and column names for clarity
  rownames(dist_matrix) <- rownames(alignment)
  colnames(dist_matrix) <- rownames(alignment)
  
  return(dist_matrix)
}

scale_tree <- function(tree, root_heigh){
  
  current_height <- max(node.depth.edgelength(tree))
  
  # 3. Compute scaling factor
  scaling_factor <- 25 / current_height
  
  # 4. Rescale branch lengths
  tree$edge.length <- tree$edge.length * scaling_factor
  
  return(tree)
}

plot_scaled_tree <- function(tree, root_height){
  
  # Plot the rescaled tree with time axis
  plot(tree,
       main = paste0("UPGMA Tree (rescaled to height ", root_height, ")"),
       x.lim = c(0, root_height),          # set x-axis range explicitly
       show.tip.label = TRUE,     # show labels on tips
       direction = "rightwards",  # default direction
       cex = 1                    # font size
  )
  
  # Add axis with time units
  axis(1, at = seq(0, 25, by = 5), labels = seq(0, 25, by = 5))
}

## function to read-in the simulated alignments
parseSimulatedTapeAlignment <- function(seed, path_to_data) {
  full_data <- data.frame()
  for(tapeNR in 1:10) {
    file_name <- paste0(path_to_data, "/simulate_alignment_and_tree.seed=",seed,".",tapeNR,".alignment.nexus")
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

build_upgma_tree = function(edit_table, nSites, collapsedSites=T){
  
  
  # edit_cell_table_65 = Ordered cell-by-65EditSites table, including non-existing sites before contraction
  if("nUMI" %in% colnames(edit_table)){
    edit_cell_table_65 <- select(edit_table, -nUMI)
  }else{
    edit_cell_table_65 <- edit_table
  }
  
  edit_cell_table_65 <- edit_cell_table_65 %>%
    pivot_longer(cols = c('Site1','Site2','Site3','Site4','Site5'), names_to = 'Sites', values_to ='Insert') %>%
    pivot_wider(id_cols = Cell, names_from = c(TargetBC,Sites), names_sep = ".", values_from = Insert)
  
  edit_cell_table_65 <- arrange(edit_cell_table_65,Cell) %>%
    select(order(colnames(edit_cell_table_65)))
  
  
  ###########
  
  sub_edit65 <- as.matrix(select(edit_cell_table_65,-Cell))
  sub_edit65[is.na(sub_edit65)] <- 'None'
  rownames(sub_edit65) <- edit_cell_table_65$Cell
  
  if(collapsedSites){
    
  }
  #sub_edit59 <- as.matrix(select(edit_cell_table_65,-c('Cell','TGGACGAC.Site5','TTTCGTGA.Site5','TGGTTTTG.Site5',
  #                                                     'TTCACGTA.Site3','TTCACGTA.Site4','TTCACGTA.Site5')))
  #rownames(sub_edit59) <- edit_cell_table_65$Cell
  cell_list <- edit_cell_table_65$Cell
  
  
  # =====================================================================
  # Generating the phylogenetic tree based on edits
  # =====================================================================
  
  shared_edit_matrix <- fun_shared_edit_matrix(sub_edit65)
  shared_edit_matrix <- as.matrix(shared_edit_matrix)
  diag(shared_edit_matrix) <- nSites
  
  distance_matrix <- nSites - shared_edit_matrix # Phylogenetic distance caludated as (# of possible sites - # of shared sites)
  distance_matrix <- as.matrix(distance_matrix)
  tree <- as.phylo(hclust(as.dist(distance_matrix), "average")) # tree built using UPGMA
  
}
# Function for calculating shared_edit_matrix
# shared_edit_matrix = Counting all shared edits per cell-pair, consistent with the sequential editing on DNA Tape
# From Choi
fun_shared_edit_matrix <- function(x) {
  #if (1 == 1){
  sub_edit65 <- x
  sub_edit65[sub_edit65 == 'None'] <- 1:filter(as.data.frame(table(sub_edit65)), sub_edit65 == 'None')$Freq
  ncell <- nrow(sub_edit65)
  cell_list <- sort(rownames(sub_edit65))
  shared_edit_matrix <- matrix(0, ncell,ncell)
  colnames(shared_edit_matrix) <- cell_list
  rownames(shared_edit_matrix) <- cell_list
  for (ii in 1:(ncell)){
    cell1 <- cell_list[ii]
    for (jj in (ii):ncell){
      cell2 <- cell_list[jj]
      for (kk in seq(0,(dim(sub_edit65)[2]-5),5)){
        if (sub_edit65[cell1,(kk+1)] == sub_edit65[cell2,(kk+1)]){
          shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
          if (sub_edit65[cell1,(kk+2)] == sub_edit65[cell2,(kk+2)]){
            shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
            if (sub_edit65[cell1,(kk+3)] == sub_edit65[cell2,(kk+3)]){
              shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
              if (sub_edit65[cell1,(kk+4)] == sub_edit65[cell2,(kk+4)]){
                shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
                if (sub_edit65[cell1,(kk+5)] == sub_edit65[cell2,(kk+5)]){
                  shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
                }
              }
            }
          }
        }
      }
      shared_edit_matrix[jj,ii] <- shared_edit_matrix[ii,jj]
    }
  }
  return(shared_edit_matrix)
}
