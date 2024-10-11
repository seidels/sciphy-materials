## ---------------------------
##
## Script purpose: Draw simulation parameters from prior
## distribution for well-calibrated analysis.
##
##
## Author: Sophie Seidel
##
## Date Created: 2023-03-31
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


## specify output location
outputDir = paste0("simulation_parameters/")

if (! dir.exists(outputDir)){
  dir.create(outputDir)
}


for (seed in 1:100){

  # 1. Substitution Model
  set.seed(seed)

  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))

  clockRate = round(rlnorm(n = 1, meanlog = -2, sdlog = 0.5), digits = 3)


  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate")

  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_", seed, ".csv"),

            quote = F, row.names = T,)
}
