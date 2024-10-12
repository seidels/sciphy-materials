#Note! tape(s) and TargetBC(s) are used interchangeably here!
get_cell_by_tape_matrix = function(dat_sanjay){

  tapes_and_edits_for_cells = data.frame(Cell = "", TargetBC = "",
                                         site1 = "", site2 = "", site3="",
                                         site4="", site5="", site6="") # Init with first dummy row

  for (cell_id in 1:nrow(dat_sanjay)){
    print(paste("=== Getting TargetBCs and edits for cell number ", cell_id, "===="))

    cell = dat_sanjay[cell_id, 1]
    cell_dat = dat_sanjay[cell_id, ]
    tapes_and_edits_for_cell = get_tapes_for_cell(cell_dat)
    tapes_and_edits_for_cells = rbind(tapes_and_edits_for_cells, tapes_and_edits_for_cell)
  }

  tapes_and_edits_for_cells = tapes_and_edits_for_cells[2: nrow(tapes_and_edits_for_cells), ] #Remove first dummy row

  return(tapes_and_edits_for_cells)
}



get_tapes_for_cell = function(cell_dat){

  non_zero_entries = cell_dat[,  which(cell_dat != "None")]
  cell_label = non_zero_entries[1, 1]
  tape_data = non_zero_entries[, non_zero_entries != cell_label]

  unique_tapes = unique(unname(sapply(names(tape_data), function(x){
    strsplit(x, split = ".Site", fixed = T)[[1]][1]
  })))

  tapes_and_edits_for_cell = data.frame(Cell = cell_label, TargetBC=unique_tapes,
                                        site1 = "", site2 = "", site3="",
                                        site4="", site5="", site6="")

  for( tape_id in 1:nrow(tapes_and_edits_for_cell)){

    #print(tape_id)
    tape = tapes_and_edits_for_cell[tape_id, "TargetBC"]
    edits = get_edits_for_tape(tape, tape_data = cell_dat)

    tapes_and_edits_for_cell[tape_id, paste0("site", 1:6)] = edits
  }

  return(tapes_and_edits_for_cell)
}

get_edits_for_tape = function(tape, tape_data){

  all_tapes_and_sites = names(tape_data)
  tape_specific_sites = startsWith(x = all_tapes_and_sites, prefix = tape)
  edits_for_tape = tape_data[, which(tape_specific_sites)]

  return(unname(unlist(edits_for_tape)))
}

