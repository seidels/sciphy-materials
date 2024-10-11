# Run well-calibrated simulations

1. Draw editing model parameters from prior distribution
2. Simulate alignments given tree and editing model parameters
3. Infer editing model parameters
4. Evaluate inference by comparing simulated (true) editing model parameters to estimated parameters

Steps 1-3 are called by the run.sh script.

Once the inferences are all converged, perform step 4 by calling
`get_summary_from_all_logs.R`

We provide exemplary files (usually for see=1) of the output in the directories simulation_parameters, simulated_data, inference_logs.