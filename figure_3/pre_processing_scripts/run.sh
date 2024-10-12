!#/bin/bash

# Filter the data
Rscript filter_data.R

# Generate alignments to be found in pre_processed_data
Rscript sequence_process_xml.R
