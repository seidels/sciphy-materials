#README

TAPE_10X_read2fromBAM_v2.py is used to get all TAPE sequencing reads from the CellRanger’s BAM output

BamExtractV2_postprocess.py is used to process them (selecting a single most frequent editing pattern per integrant in each cell and getting the edits by looking at GGA key sequences on each sites).

To generate the entropy ranked list of TargetBCs per individual gastruloid the script ‘mGASv2_TargetBC_filtering_v2.R’ is used.