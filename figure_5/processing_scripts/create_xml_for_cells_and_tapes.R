## ---------------------------
##
## Purpose of script: Create xml for selected cells and tapes
##
## Author: Sophie Seidel
##
## Date Created: 2022-12-13
##
## Copyright (c) Sophie Seidel, 2022
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------
# load libs

library(stringr)
library(dplyr)

## ---------------------------
source("processing_scripts/functions.R")

## ---------------------------

output_folder = "processed_data/"

#load the filtered data:

edit_table = readRDS("processed_data/cell_by_tape.RDS")

site_columns = 3:8

# 2nd, create a map from the trinucleotides to an integer
insert_to_integer_map = get_insert_to_integer_map(edit_table, columns=site_columns, dataset = "gastruloid")
write.csv(insert_to_integer_map, file = "processed_data/1-insert_to_integer_map.csv")

# convert these trinucleotides to integer in the edit table
for(site_column in site_columns) {
  edit_table[, site_column] <- unlist(lapply(edit_table[, site_column], function(x) {pick_code(x, insert_to_integer_map)}))
}

# concatenate edits at sites, adding commas:
edit_table$beast_seq <- apply(edit_table[, site_columns], 1, function(x) {str_flatten(x,collapse = ",")})


# Here for the first set of filtered cells and tapes
targetBCs = readRDS("processed_data/1-15toptapes.RDS")
cells = readRDS("processed_data/1-cells_15toptapes_75cells.RDS")

edit_table_filtered = edit_table[(edit_table$Cell %in% cells) & (edit_table$TargetBC %in% targetBCs), ]
write.csv(x = edit_table_filtered, file = paste0(output_folder, "1-edit_table.csv"))
write_alignment_to_xml(cell_sample = cells, dataset = edit_table_filtered, targetBCs = targetBCs, n_cells = 75,
                           alignment_file = paste0(output_folder, "1-alignment.txt"), n_states = 39)



