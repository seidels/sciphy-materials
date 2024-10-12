#!/bin/bash

# Run the data pre-processing pipeline

# First take input matrix and create our desired input matrix format,
# namely a cell by tape matrix
Rscript create_cell_by_tape_matrix.R

# Second, filter the data
Rscript select_cells.R

# Third, create the alignment in xml format
Rscript create_xml_for_cells_and_tapes.R

