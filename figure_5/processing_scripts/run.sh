#!/bin/bash

# Run the data pre-processing pipeline

# First, filter the data by 8 pre-selected barcodes and keep only cells that
# have all of them
Rscript filter_by_selected_barcodes.R


# Second, create the alignment in xml format
Rscript write_alignments.R

