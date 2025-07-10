

# General description

Figure 4 compares estimates using SciPhy with ones obtained with UPGMA on the HEK293T dataset and shows a benchmark of SciPhy against existing phylogenetic tools.

## SciPhy UPGMA comparison on real data

1) Panel A + Appendix Figures 3-6: Visualising the differences in tree estimates across both methods by mapping distances between UPGMA and SciPhy estimates in 2D, using R packages TreeDist and treespace. 
2) Panel B: Visualising differences in the corresponding branch lengths by plotting the LTT plot for the estimated trees.
3) Panel C: Visualising the downstream effect on the growth rates that may be estimated using such trees.

All 3 panels are visualized by running `plot_figure_4.R`, after generating input files as described below.

# Panel A 

To generate input files for the plots in panel A:

First, generate the inputs for TreeDist and treespace:

1) Thin down the tree traces obtained from running `figure_3/inference_output/3-typewriter_sampling_clockPerSite_13Sites_1000Cells_DataSet1.xml`. To do so, run: `logcombiner -resample 3000000 -log combined_clockPerTarget_sampling_DataSet1.log -o clockPerTarget_sampling_DataSet1_3000000.log` and `logcombiner -resample 3000000 -log combined_clockPerTarget_sampling_DataSet1.log -o clockPerTarget_sampling_DataSet1_3000000.trees` and place these files in `inference_output` This is needed as mapping distances between the (very large) full posterior set of trees is not tractable.
2) Build a UPGMA tree on the alignment. For this, run `build_UPGMA_Dataset1.R`. 

Second, calculate distances between these sets of trees and map these distances in 2D by running : `run_TreeDist_CI.R`. (Also `run_TreeDist_PI.R`, `run_TreeDist_RF.R` to generate outputs for the supplemental figures). 

# Panel B

Panel B relies on: 

1) Running a SciPhy analysis on a fixed UPGMA topology to estimate its branch lengths (along with the growth rate) the corresponding xml file to run is: `inference_output/1000_UPGMA_medianPosteriorHeight_estimateBranchLengths_infer_rho_sampling.xml`. 
2) Generating and plotting LTT plots for the UPGMA trees, the SciPhy estimated trees, and the UPGMA with SciPhy estimated branch lengths.

# Panel C

To generate panel C, 2 more SciPhy analyses that estimate the growth rate on a fixed UPGMA topology are needed: 
1) `inference_output/1000_UPGMA_infer_rho_sampling.xml`
2) `inference_output/1000_UPGMA_infer_rho_sampling_median_height.xml`

# Visualising results

Once all of these scripts have been run and their outputs are in `inference_output`, the entire figure 4 is reproduced by running `plot_figure_4.R`. 
Appendix figures 3-6 are generated with the outputs described for Panel A by running `create_supp_figs_3_4_5.R` and `create_supp_figs_6.R` in supplemental_figures. 

## SciPhy vs existing methods benchmark
All under  benchmark_simulated_data.