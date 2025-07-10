mark_unedited_sites = function(edit_table){
  # set unedited sites to 'None' , s.t. this is treated as a shared site in upgma
  edit_table[is.na(edit_table)] <- 'None'
  
  return(edit_table)
}

set_non_contracted_sites_to_na = function(edit_table){
  
  edit_table[edit_table$TargetBC == 'TGGACGAC',7] <- NA
  edit_table[edit_table$TargetBC == 'TTTCGTGA',7] <- NA
  edit_table[edit_table$TargetBC == 'TGGTTTTG',7] <- NA
  edit_table[edit_table$TargetBC == 'TTCACGTA',5:7] <- NA
  
  return(edit_table)
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