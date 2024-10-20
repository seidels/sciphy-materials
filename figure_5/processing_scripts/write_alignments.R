## ---------------------------
##
## Script name: write alignments
##
## Purpose of script: write alignments given targets per cell input data
##
## Author: Sophie Seidel
##
## Date Created: 2023-04-19
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

library(plyr)

## load up our functions into memory

source("useful_scripts_across_categories.R")

## ---------------------------

# output files
alignment_file = "../processed_data/alignment_filtered_for_8barcodes.xml"
integers_dat_file = "../processed_data/edits_to_integer_map.csv"
integers_conversion_table = "../processed_data/integer_conversion_table.csv"
traits_file = "../processed_data/date_traits.xml"

# input file
filtered_dat_file = "../processed_data/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS"
filtered_dat = readRDS(filtered_dat_file)

# pre-process
conversion <- convert_edits_to_integer_with_edit_list(filtered_dat, number_of_sites = 5)
# get the converted sequences
integers_dat <- conversion$targets_per_cell
# get the mapping of trinucleotide to integer
conversion_table <- conversion$possible_edits

# save this mapping:
write.csv(x = conversion_table, file = integers_conversion_table)

## Remove because Site 6 is unedited throughout
integers_dat = integers_dat[, ! names(integers_dat) %in% c("Site6")]
write.csv(x = integers_dat, file = integers_dat_file)

# write alignment
for (targetBC in unique(integers_dat$TargetBC)){

  print(targetBC)
  write_targetBC_alignment_to_xml(targets_per_cell_dat = integers_dat,
                                  filename = alignment_file,
                                  targetBC = targetBC)
}

# date trait
taxon_and_date_individual = unname(sapply(filtered_dat$Cell, function(x){
  paste(x, "=10", sep = "")
}))
# to be pasted in xml
taxon_and_date_merged = paste(taxon_and_date_individual, collapse = ",")
write(taxon_and_date_merged, file = traits_file)
