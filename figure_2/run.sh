!#/bin/bash

## This script runs the validation study

# requirements:
## BEAST v2.7, sciphy v1.0.0 and bdmm prime installed
## Change all_beasts_2.7_2024.jar (see below) to your BEAST executable where sciphy is also installed
## R version 4.4.1 and extraDistr_1.10.0 installed 


# 1. Draw 100 editing model parameter vectors, i.e., clock rate and the insert probabilities from prior distributions
Rscript draw_simulation_params.R

# 2. Simulate 100 trees and alignments
mkdir simulated_data
for seed in `seq 1 100`
do
    java -jar ~/all_beasts_2.7_2024.jar -seed $seed -version_file version.xml simulate_alignment_and_tree.xml
done


# 3. Run 100 inferences
mkdir inference_logs
for seed in `seq 1 100`
do
    java -jar ~/all_beasts_2.7_2024.jar -seed $seed -version_file version.xml infer_given_fixed_tree_13_inserts.xml
done

    
## When running this script on the cluster, consider submitting a jobarray such as in the examples below

#sbatch --job-name="tree_alignment_simulations" --array=1-100 --wrap="java -jar ~/all_beasts_2.7_2024.jar -seed \$SLURM_ARRAY_TASK_ID -version_file version.xml simulate_alignment_and_tree.xml"

#run inference
#sbatch --job-name="validations" --array=1-100 --wrap="java -jar ~/all_beasts_2.7_2024.jar -seed \$SLURM_ARRAY_TASK_ID -version_file version.xml infer_given_fixed_tree_13_inserts.xml"
