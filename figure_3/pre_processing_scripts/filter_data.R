## ---------------------------
##
## Script name: filtering.R
##
## Purpose of script: Filter out cells based on the following criteria.
## - they do not have all of the 13 most frequently occurring TargetBCs
##
## Author: Sophie Seidel
##
## Date Created: 2022-08-22
##
## Copyright (c) Sophie Seidel, 2022
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes: data obtained from https://github.com/shendurelab/DNATickerTape/
##
##
## ---------------------------

# load the preprocessed data ---------

edit_table = read.csv("../data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F,
                      header = T, na.strings=c("","NA"))

#length(unique(edit_table$Cell)) == 16823

## filtering
### filter out entries that have no TargetBC or UMI
edit_table = edit_table[which(!(is.na(edit_table$TargetBC))), ]
edit_table = edit_table[which(!(is.na(edit_table$nUMI))), ]
edit_table = edit_table[which(!(is.na(edit_table$Cell))), ]
#length(unique(edit_table$Cell)) == 16810


### filter out entries that have not correctly ordered edits
#### pairwise comparisons, remove all entries without edit at site 1 and edit at site 2-5
edit_table = edit_table[! (is.na(edit_table[, 3]) &
                             !(is.na(edit_table[, 4]))), ]
edit_table = edit_table[! (is.na(edit_table[, 3]) &
                             !(is.na(edit_table[, 5]))), ]
edit_table = edit_table[! (is.na(edit_table[, 3]) &
                             !(is.na(edit_table[, 6]))), ]
edit_table = edit_table[! (is.na(edit_table[, 3]) &
                             !(is.na(edit_table[, 7]))), ]
#### same for site 2
edit_table = edit_table[! (is.na(edit_table[, 4]) &
                             !(is.na(edit_table[, 5]))), ]
edit_table = edit_table[! (is.na(edit_table[, 4]) &
                             !(is.na(edit_table[, 6]))), ]
edit_table = edit_table[! (is.na(edit_table[, 4]) &
                             !(is.na(edit_table[, 7]))), ]

#### same for site 3
edit_table = edit_table[! (is.na(edit_table[, 5]) &
                             !(is.na(edit_table[, 6]))), ]
edit_table = edit_table[! (is.na(edit_table[, 5]) &
                             !(is.na(edit_table[, 7]))), ]

#### same for site 4
edit_table = edit_table[! (is.na(edit_table[, 6]) &
                             !(is.na(edit_table[, 7]))), ]

# length(unique(edit_table$Cell)) == 16810

### get the most frequent TargetBCs && assert that they are the same found in the original paper
frequency_of_target_bcs = as.data.frame(table(edit_table$TargetBC, useNA = "ifany"))
frequency_of_target_bcs = frequency_of_target_bcs[order(frequency_of_target_bcs$Freq, decreasing = T), ]
frequent_target_bcs = frequency_of_target_bcs[1:13, "Var1"]

targetBCs_from_paper = c("ATGGTAAG", "ATTTATAT",
                         "ATTTGGTT", "GCAGGGTG",
                         "GTAAAGAT", "TAGATTTT",
                         "TGCGATTT", "TGGACGAC",
                         "TGGTTTTG", "TTAGATTG",
                         "TTCACGTA", "TTGAGGTG",
                         "TTTCGTGA")

#sort(frequent_target_bcs)
# Make sure that we get the same 13 TargetBCs as Choi in his paper
assertthat::are_equal(length(intersect(frequent_target_bcs, targetBCs_from_paper)), 13)


### filter out cells that do not have the 13 most frequent TargetBCs && assert that no cell has the same targetBC twice
for (cell in unique(edit_table$Cell)){

  print(cell)

  cell_target_bcs = edit_table[which(edit_table$Cell == cell), "TargetBC"]

  if (length(cell_target_bcs) != length(unique(cell_target_bcs))){
    stop("Same target bc in one cell")
  }

  keep_cell = all(frequent_target_bcs %in% cell_target_bcs)

  if (! keep_cell){
    edit_table = edit_table[which(!(edit_table$Cell == cell)), ]
  }
}
## Here, we have 3221 cells which is close to the reported 3257
#length(unique(edit_table$Cell)) == 3221

# subset edit table to those targetBCs that we will continue to use
## as expected nrow(edit_table) = 13 (targetBCs) * 3221 (cells)
edit_table = edit_table[edit_table$TargetBC %in% frequent_target_bcs, ]

## Here, we have 3221 cells
saveRDS(object = edit_table, file = "../data/edit_table_filtered.RDS")