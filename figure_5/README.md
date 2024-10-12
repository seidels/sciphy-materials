# Figure 5

To generate figure 5C, the filtered data can be plotted using the script 


To generate the plots B, D-E in figure 5, the analysis on 2 time change points detailed in `inference_output/3-mGASv2-skyline-ou_1000000.xml` has to be run.

Then run TODO , to generate figures B and E, and `plot_figure_5_D.R` to generate figure 5D.

## Generate alignments

**Please note that the data will only be made publicly available upon publication.** In the meantime, we are sharing empty placeholder files, such as `processed_data/cell_by_tape.RDS`, to demonstrate the structure of the data pipeline.


To generate the alignment underlying the analyses, run the `run.sh` executing the data preprocessing pipeline to be found in  `./processing_scripts/`.

## Supplemental figures
To generate Appendix 0—figure 12., the analysis with 3 time change points has to be run, i.e. `inference_output/4-mGASv2-skyline-ou_1000000.xml

Then run `plot_supp_fig_12.R` to generate supplemental figure 12.
