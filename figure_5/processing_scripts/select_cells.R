
cell_by_tape = readRDS(file = "processed_data/cell_by_tape.RDS")


library(ggplot2)

# Filtering: We hope to find an alignment with at least x tapes per cell
# 1) Thus, filter out all cells that do not have at least x tapes
# 2) and filter out all tapes that are not present in more than x cells

x = 40

## ---------------------------
## 1) Filter cells
cells = unique(cell_by_tape$Cell)
cell_by_n_tapes = data.frame(cell = cells, n_tapes = 0)

for (cell_id in 1:length(cells)){
  cell = cells[cell_id]
  n_tapes = nrow(cell_by_tape[cell_by_tape$Cell == cell, ])
  cell_by_n_tapes[cell_id, "n_tapes"] = n_tapes
}

cell_by_n_tapes = cell_by_n_tapes[order(cell_by_n_tapes$n_tapes), ]

## filtering out cells that have less than x tapes
cells_with_too_few_tapes = cell_by_n_tapes[which(cell_by_n_tapes$n_tapes < x), "cell"]

## ---------------------------
## 1) Filter barcodes
tapes = unique(cell_by_tape$TargetBC)
tape_by_n_cells = data.frame(tape=tapes, n_cells = 0)

for (tape_id in 1:nrow(tape_by_n_cells)){
  tape = tapes[tape_id]
  n_cells = nrow(cell_by_tape[cell_by_tape$TargetBC == tape, ])
  tape_by_n_cells[tape_id, "n_cells"] = n_cells
}

tape_by_n_cells = tape_by_n_cells[order(tape_by_n_cells$n_cells), ]
ggplot(tape_by_n_cells, aes(x=n_cells))+
  geom_histogram(bins = 10)

tapes_with_too_few_cells = tape_by_n_cells[which(tape_by_n_cells$n_cells < x), "tape"]

dat_filtered = cell_by_tape[which(!(cell_by_tape$Cell %in% cells_with_too_few_tapes)), ]
dat_filtered = dat_filtered[which(!(cell_by_tape$TargetBC %in% tapes_with_too_few_cells)), ]

## ---------------------------
length(unique(dat_filtered$Cell))
length(unique(dat_filtered$TargetBC))

## Problem - filtering did not reduce the number such that exhaustive search is feasible.

## ---------------------------
# Take the top y tapes and tell me in how many cells they all are
y = 10
tape_by_n_cells = tape_by_n_cells[order(tape_by_n_cells$n_cells, decreasing = T), ]
tapes_y = tape_by_n_cells[1:y, "tape"]

cells = unique(cell_by_tape$Cell)

cells_with_y_tapes = c()

for (cell_id in 1:length(cells)){

  cell = cells[cell_id]
  tapes_in_cell = cell_by_tape[which(cell_by_tape$Cell == cell), "TargetBC"]

  if (all(tapes_y %in% tapes_in_cell)){
    cells_with_y_tapes = c(cells_with_y_tapes, cell)
  }
}
cells_with_y_tapes
saveRDS(object = cells_with_y_tapes, file = "processed_data/1-cells_15toptapes_75cells.RDS")
saveRDS(object = tapes_y, file = "processed_data/1-15toptapes.RDS")

