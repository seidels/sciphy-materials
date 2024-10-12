## ---------------------------
##
## Purpose of script: Take the input matrix provided by sanjay and convert it into a
## cell by tape matrix
##
## Author: Sophie Seidel
##
## Date Created: 2024-04-12
##
## Copyright (c) Sophie Seidel, 2024
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## load up our functions into memory

source("processing_scripts/functions.R")

## ---------------------------

dat = "data/input_matrix.csv"

dat = read.csv(file = dat)

### sanjay dat ###

### rows are cells and the rownames give the cell barcode

### columns are individual Tape sites that were edited.
### For example, if the column name is "AAAGTAAATCTC-0-Site1", it means
### AAAGTAAATCTC was the Tape barcode, 0 is the specific locus of integration
### for that barcode, and we're looking at the first of the 6 possible sites
### that could have been edited for that locus.


cell_by_tape_dat = get_cell_by_tape_matrix(dat)

saveRDS(cell_by_tape_dat, file = "cell_by_tape.RDS")

