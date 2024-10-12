#!/bin/bash

#load java modules
env2lmod
module load openjdk/17.0.0_35

sbatch --time=120:00:00 --job-name="cell_culture_2_sites" --array=1-3 --wrap="java -jar ~/all_beasts_2.7_2024.jar -threads 6 -seed \$SLURM_ARRAY_TASK_ID -version_file version.xml analysis_2sites_DataSet1.xml"
