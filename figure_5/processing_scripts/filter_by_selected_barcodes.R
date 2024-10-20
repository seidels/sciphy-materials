## ---------------------------
##
## Script name: Filter by selected barcodes
##
## Purpose of script: Keep only cells that have all selected
## barcodes.
##
## Author: Sophie Seidel
##
## Date Created: 2023-04-18
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

## load up our functions into memory

source("useful_scripts_across_categories.R")

## ---------------------------

output_file = "../processed_data/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS"

## read files

dat_file = "../data/mGASv2_Lane2_CellByTape_10X_bamExtractV2_t3_collapse.csv"
dat = read.csv(dat_file)

selected_barcode_file = "../data/mGASv2_TargetBC_selected8.csv"
selected_barcodes_dat = read.csv(selected_barcode_file)
selected_barcodes = selected_barcodes_dat$TargetBC


## filter to only keep the cells in dat that have all 8 selected barcodes
cells_in_dat = unique(dat$Cell)
cells_filtered =  sapply(cells_in_dat, FUN = function(x){
  cell_has_selected_barcodes(cell = x, dat = dat,
                             selected_barcodes = selected_barcodes)

})

cells_filtered = cells_in_dat[cells_filtered]
dat_filtered_cells = dat[which(dat$Cell %in% cells_filtered), ]

dat_filtered_cells_barcodes = dat_filtered_cells[which(dat_filtered_cells$TargetBC %in% selected_barcodes), ]


## save output
saveRDS(object = dat_filtered_cells_barcodes, file = output_file)

