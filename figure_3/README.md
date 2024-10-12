

# SciPhy analyses

Figure 3 builds on 3 analyses.

1) Estimating one growth rate over the experiment, while the sampling propotion is set to the experimentally reported one.
`1-typewriter_clockPerSite_13Sites_1000Cells_DataSet1.xml`

2) Estimating a dynamic (or time-varying) growth rate over the experiment, while the sampling proportion is set to the experimentally reported one.
`2-typewriter_SKY_OU_clockPerSite_13Sites_1000Cells_DataSet1.xml`

3) Estimating one growth rate over the experiment, while co-estimating the sampling proportion. TODO

Each analysis xml can be found under the `inference_output` directory. 

# Generate alignment

The 3 analyses use the same underlying data. To filter the data, located in `data`, pre-process it and generate the alignments, run the `pre_processing_scripts/run.sh` script.

You can find the alignment under `pre_processed_data/1-alignment_seed1.txt`. 

# Visualising results

Figure 3 B-D and F-G build on analysis (1). 
To reproduce figure 3 B-D, run the `plot_inference_results.R` script and for F-G use the `plot_mcc_tree_annotated.R` scripts.

To reproduce figure 3, E, which builds on analysis (3), run TODO.
