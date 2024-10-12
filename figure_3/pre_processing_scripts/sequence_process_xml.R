## ---------------------------
##
## Script name: sequence_process_xml
##
## Purpose of script: Pre-processing and sampling of sequences for xml inputs
##
## Author: Antoine Zwaans & Sophie Seidel
##
## Date Created: 2022-12-13
##
## Copyright (c) Antoine Zwaans & Sophie Seidel, 2022
## Email: antoine.zwaans@bsse.ethz.ch, sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## ---------------------------
# load libs

library(stringr)
library(dplyr)

## ---------------------------
source("functions.R")

## ---------------------------

output_folder = "../pre_processed_data/"

#load the filtered data:

edit_table_by_5 = readRDS("../data/edit_table_filtered.RDS")

site_columns = 3:7

# 1st Convert 6-mers to 3-mers
# by removing the last 3 nucleotides that do not contain information
# (they are always GGA); replace <Na> by "NA"

for(site_column in site_columns) {

  edit_table_by_5[, site_column] <- substring(edit_table_by_5[, site_column], 1,3)
  edit_table_by_5[which(is.na(edit_table_by_5[, site_column])), site_column] <- "NA"

}

# 2nd, create a map from the trinucleotides to an integer
insert_to_integer_map = get_insert_to_integer_map(edit_table_by_5, columns=site_columns, dataset = "cell_culture")
n_states = length(insert_to_integer_map$integer)
write.csv(insert_to_integer_map, file = "../pre_processed_data/1-insert_to_integer_map.csv")

# convert these trinucleotides to integer in the edit table
for(site_column in site_columns) {
  edit_table_by_5[, site_column] <- unlist(lapply(edit_table_by_5[, site_column], function(x) {pick_code(x, insert_to_integer_map)}))
}

# concatenate edits at sites, adding commas:
edit_table_by_5$beast_seq <- apply(edit_table_by_5[, site_columns], 1, function(x) {str_flatten(x,collapse = ",")})

#targetBCs as used in the paper:
targetBCs = c("ATGGTAAG","ATTTATAT",
                         "ATTTGGTT", "GCAGGGTG",
                         "GTAAAGAT", "TAGATTTT",
                         "TGCGATTT", "TGGACGAC",
                         "TGGTTTTG", "TTAGATTG",
                         "TTGAGGTG",
                         "TTTCGTGA","TTCACGTA")

# create alignments for different cell numbers, 3 different sub samples each
n_cells = 1000
seed = 1

output_folder_ncells = paste0(output_folder)

cell_sample = subsample_dataset(n_cells, edit_table_by_5, seed)

edit_table_for_sample = edit_table_by_5[ (edit_table_by_5$Cell %in% cell_sample) & (edit_table_by_5$TargetBC %in% targetBCs), ]
write.csv(x = edit_table_for_sample, file = paste0(output_folder_ncells, "1-edit_table_sample_", seed, ".csv"))

write_cell_ids_to_file(cell_sample = cell_sample, cell_ids_file = paste0(output_folder_ncells, "1-cell_ids_seed", seed, ".txt"))
write_alignment_to_xml(cell_sample = cell_sample, dataset = edit_table_by_5, targetBCs = targetBCs, n_cells = n_cells,n_states = n_states,
                           alignment_file = paste0(output_folder_ncells, "1-alignment_seed", seed, ".txt"))

print(paste("Successfully generated alignments for 1000 cells"))


