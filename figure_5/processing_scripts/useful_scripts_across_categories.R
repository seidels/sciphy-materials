
library(plyr)

get_nr_of_edits = function(edit_table, unedited_symbol=NA){

  if(!(is.na(unedited_symbol))){
    edit_table[edit_table == unedited_symbol] = NA
  }
  ## compare to data
  nr_of_edits = rep(x = 0, nrow(edit_table))

  for (i in 1:nrow(edit_table)){

    if(is.na(edit_table[i, "Site1"])){

      nr_of_edits[i] = 0
      next

    }else if(is.na(edit_table[i, "Site2"])){

      nr_of_edits[i] = 1
      next

    }else if (is.na(edit_table[i, "Site3"])){

      nr_of_edits[i] = 2
      next

    }else if (is.na(edit_table[i, "Site4"])){

      nr_of_edits[i] = 3
      next

    }else if (is.na(edit_table[i, "Site5"])){

      nr_of_edits[i] = 4
      next

    }else{
      nr_of_edits[i] = 5
    }
  }
  return(nr_of_edits)
}


cell_has_selected_barcodes = function(cell, dat, selected_barcodes){

  cell_barcodes = dat[dat$Cell == cell, "TargetBC"]

  if (all(selected_barcodes %in% cell_barcodes)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

convert_edits_to_integer = function(targets_per_cell, number_of_sites=5){

  possible_edits = get_possible_edits(targets_per_cell)

  for (site_i in 1:number_of_sites){

    site_string = paste0("Site", site_i)
    targets_per_cell[, site_string] = mapvalues(x = targets_per_cell[, site_string], from = possible_edits$edits, to = possible_edits$integer)
  }

  return(targets_per_cell)
}

convert_edits_to_integer_with_edit_list = function(targets_per_cell, number_of_sites=5){
  
  possible_edits = get_possible_edits(targets_per_cell)
  
  
  
  for (site_i in 1:number_of_sites){
    
    site_string = paste0("Site", site_i)
    targets_per_cell[, site_string] = mapvalues(x = targets_per_cell[, site_string], from = possible_edits$edits, to = possible_edits$integer)
  }
  
  return(list(possible_edits=possible_edits,targets_per_cell=targets_per_cell))
}

get_possible_edits = function(targets_per_cell){

  columns_with_edits = colnames(targets_per_cell)
  columns_with_edits = columns_with_edits[startsWith(x = columns_with_edits, prefix = "Site")]

  all_edits = unique(unlist(targets_per_cell[ , columns_with_edits]))
  all_edits = all_edits[all_edits != "None"]

  edit_to_integer_map = data.frame(edits = all_edits, integer=1:length(all_edits))
  edit_to_integer_map[nrow(edit_to_integer_map) + 1, ] = c("None", 0)

  return(edit_to_integer_map)
}

write_targetBC_alignment_to_xml = function(targets_per_cell_dat, targetBC, filename, barcode, datName){

  write_alignment_start(targetBC, filename)

  write_sequence_into_alignment(targets_per_cell_dat, targetBC, filename)

  write_alignment_closure(filename)
}

write_alignment_start = function(targetBC, filename){

  write(x = c(
    paste0('<data  id="preliminary_gastruloid_targetBC_', targetBC,
           '" spec="Alignment" name="alignment" dataType="integer">')),

    file = filename, ncolumns = 1, append = TRUE, sep = " ")
}


write_sequence_into_alignment = function(targets_per_cell_dat, targetBC, filename){

  dat_with_targetBC = targets_per_cell_dat[which(targets_per_cell_dat$TargetBC == targetBC), ]

  cells = unique(dat_with_targetBC$Cell)

  for (i in 1:length(cells)){

    cell  =  cells[i]
    targets_per_cell = dat_with_targetBC[which(dat_with_targetBC$Cell==cell), ]

    targets_concatenated = paste(
      targets_per_cell[startsWith(colnames(targets_per_cell), prefix = "Site")], collapse = ","
    )

    write(x = c(paste0('<sequence id="', (i-1), '_cell_', cell, '_targetBC_', targetBC,
                       '" spec="Sequence" taxon="',
                       cell,'" value="', targets_concatenated, '"/>')),
          file = filename, ncolumns = 1, append = TRUE, sep = " ")
  }
}

write_alignment_closure = function(filename){
  write(x = "</data>", file = filename, append = TRUE, sep = " ")
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
