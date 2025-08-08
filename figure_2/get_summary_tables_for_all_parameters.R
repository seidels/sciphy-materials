## ---------------------------
##
## Script name: get summary for all inferred parameters in the full validation
##
## Purpose of script: Analyse well-calibrated simulations
##
## Author: Antoine Zwaans, Sophie Seidel
##
## Date Created: 2024-07-05
##
## Copyright (c) Antoine Zwaans, Sophie Seidel, 2024
## Email: azwaans@ethz.ch, ophie.seidel@posteo.de
##
##
## ---------------------------


#load packages
library(tracerer)
library(HDInterval)
library(ggplot2)
library(beastio)
library(stringr)
library(tidyverse)
library(phytools)
library(treebalance)


log_dir = "../inference_logs"

simulation_dir = "../simulation_parameters/"

inc = function(x){
  eval.parent(substitute(x <- x+1))
}

is_parameter_in_hpd = function(hpd_lower, hpd_upper, true_parameter){
  if(true_parameter >= hpd_lower && true_parameter <= hpd_upper)
    return(TRUE)
  else
    return(FALSE)
}

## extract tree height and tree length from true simulated trees
num_tips <- c()
tree_heights <- c()
tree_lengths <- c()
tree_B1 <- c()
for(i in 1:100) {
  file_name <- paste0("simulated_data/simulate_alignment_and_tree.",i,".newick")
  tree <- readLines(file_name)
  tree <- str_remove_all(tree,"\\[.........\\]")
  tree <- str_remove_all(tree,"\\[..........\\]")
  tree <- str_remove_all(tree,"\\[...........\\]")
  tree <- paste0(tree,";")
  p_tree <- ape::read.tree(text = tree)
  tree_heights <- c(tree_heights,max(nodeHeights(p_tree)))
  tree_lengths <- c(tree_lengths,sum(p_tree$edge.length))
  tree_B1 <- c(tree_B1,treebalance::B1I(as.phylo(p_tree)))
  num_tips <- c(num_tips,length(p_tree$tip.label))
  
}

true_tree_stats <- data.frame(seed=1:100,ntips=num_tips,true_heights=tree_heights,true_lengths=tree_lengths,true_b1=tree_B1)

# extract inference results from all datasets
nr_converged_chains = 0

clock_rate_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

insert_rate_inference = data.frame(seed=rep(1:100, each = 13), hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,
                                   insertRate= rep(1:13,100))

tree_height_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

tree_length_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

for (seed in  1:100){

  print(seed)
  # get inference log
  log_file = paste0("../inference_logs/infer_given_unfixed_tree_13_inserts.", seed, ".log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.1)

  # get true parameter from simulation
  simulation_parameter_file = paste0(simulation_dir, "simParams_", seed, ".csv")
  simulation_parameters = read.csv(simulation_parameter_file)
  colnames(simulation_parameters) = c("parameter", "value")

  #check that chains has converged
  esses = calc_esses(log_data_wo_burnin, sample_interval = 1000)
  esses = esses[!colnames(esses) %in% c("birthRate", "deathRate", "prior","samplingProportion")]
  esses

  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    next

  }else{
    inc(nr_converged_chains)
  }

  ####
  # CLOCK RATE
  ###
  
  ## extract clock rate hpd
  clock_rate_hpd = hdi(object = log_data_wo_burnin$clockRate, credMass = 0.95)

  ### extract simulation value
  true_clock_rate = simulation_parameters[14, 2]

  ## check recovery
  recovered = is_parameter_in_hpd(hpd_lower = clock_rate_hpd["lower"], hpd_upper = clock_rate_hpd["upper"],
                      true_parameter = true_clock_rate)
  
  ## add all stats to table
  clock_rate_inference[seed, ] = c(seed, clock_rate_hpd, median(log_data_wo_burnin$clockRate), true_clock_rate, recovered)

  ####
  # TREE HEIGHT
  ###

  ### extrac tree height hpd
  tree_height_hpd = hdi(object = log_data_wo_burnin$'treeHeight.t.alignment', credMass = 0.95)

  ### extract simulation value
  true_tree_height =   true_tree_length = true_tree_stats[true_tree_stats$seed == seed, "true_heights"]

  ## check recovery
  recovered = is_parameter_in_hpd(hpd_lower = tree_height_hpd["lower"], hpd_upper = tree_height_hpd["upper"],
                                  true_parameter = true_tree_height)
  
  ### add all stats to table
  tree_height_inference[seed, ] = c(seed, tree_height_hpd, median(log_data_wo_burnin$'treeHeight.t.alignment'), true_tree_height, recovered)

  ####
  # TREE LENGTH
  ###
  
  ## extract tree length hpd
  tree_length_hpd = hdi(object = log_data_wo_burnin$'treeLength.t.alignment', credMass = 0.95)

  ### extract simulation value
  true_tree_length = true_tree_stats[true_tree_stats$seed == seed, "true_lengths"]

  ## check recovery
  recovered = is_parameter_in_hpd(hpd_lower = tree_length_hpd["lower"], hpd_upper = tree_length_hpd["upper"],
                                  true_parameter = true_tree_length)
  
  ### add all stats to table
  tree_length_inference[seed, ] = c(seed, tree_length_hpd, median(log_data_wo_burnin$'treeLength.t.alignment'), true_tree_length, recovered)


  ####
  # INSERT PROBABILITIES 
  ###

  for (insert_rate_nr in 1:13){

    ## extract edit probability hpd
    insert_rate = paste0("insertRates.", insert_rate_nr)
    insert_rate_hpd = hdi(object = log_data_wo_burnin[, insert_rate], credMass = 0.95)
    median_insert_rate = median(log_data_wo_burnin[, insert_rate])
    
    ## extract simulation value
    true_insert_rate = simulation_parameters[insert_rate_nr, "value"]
    
    ## check recovery
    recovered = is_parameter_in_hpd(hpd_lower = insert_rate_hpd["lower"], hpd_upper = insert_rate_hpd["upper"],
                                    true_parameter = true_insert_rate)


    row_index = which(insert_rate_inference$seed == seed & insert_rate_inference$insertRate == insert_rate_nr)
    
    ### add all stats to table
    insert_rate_inference[row_index, ] = c(seed, as.numeric(insert_rate_hpd), median_insert_rate, true_insert_rate,
                                           recovered, insert_rate_nr)
  }
}

#reorder tables

### clock rate
clock_rate_inference = clock_rate_inference[order(clock_rate_inference$true_value), ]
clock_rate_inference$orderedSeed = 1:100

### tree height
tree_height_inference = tree_height_inference[order(tree_height_inference$true_value), ]
tree_height_inference$orderedSeed = 1:100

### miss prob of alignment with length 5 targets
tree_length_inference = tree_length_inference[order(tree_length_inference$true_value), ]
tree_length_inference$orderedSeed = 1:100

#coverages per insert rate
coverages_per_insert <- c()
for(i in 1:13) {
  coverage <- sum(insert_rate_inference[which(insert_rate_inference$insertRate == i),"recovered"])/nr_converged_chains
  coverages_per_insert <- c(coverages_per_insert,coverage)
}
coverages_per_insert

#coverage for the clock rate
sum(clock_rate_inference$recovered)/nr_converged_chains
#coverage for the tree height
sum(tree_height_inference$recovered)/nr_converged_chains
#coverage for the tree length
sum(tree_length_inference$recovered)/nr_converged_chains

write.csv(x = clock_rate_inference, "inference_logs/summary/clock_rate_inference.csv",quote = F, row.names = F)
write.csv(x = tree_height_inference, "inference_logs/summary/tree_height_inference.csv",quote = F, row.names = F)
write.csv(x = tree_length_inference, "inference_logs/summary/tree_length_inference.csv",quote = F, row.names = F)
write.csv(x = insert_rate_inference, "inference_logs/summary/insert_rate_inference.csv",quote = F, row.names = F)


####
# TREE IMBALANCE
###

nr_converged_chains = 0
B1_index_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

for (seed in  1:100){
  
  print(seed)
  # get inference log
  log_file = paste0("inference_logs/infer_given_unfixed_tree_13_inserts.tree.", seed, ".treestreestats.log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.1)
  names(log_data_wo_burnin) <- c("Sample","Colless.tree.imbalance","B1")
  
  #check that chains has converged + account for different logging frequencies for TreeStat2
  esses = calc_esses(log_data_wo_burnin, sample_interval = 1000)
  if(seed %in% c(100,18,22,26,3,48,59,67,75,79,87,90,98)) {
    esses = calc_esses(log_data_wo_burnin, sample_interval = 10000)
  }
  if(seed == 2) {
    esses = calc_esses(log_data_wo_burnin, sample_interval = 20000)
  }
  
  #some statsdefault to 1, make them high ESS
  if(is.na(esses$Colless.tree.imbalance) || is.na(esses$B1)) {
    if(sum(log_data_wo_burnin$Colless.tree.imbalance) == length(log_data_wo_burnin$Colless.tree.imbalance)) {
      esses$Colless.tree.imbalance <- 400 
    }
    if(sum(log_data_wo_burnin$B1) == length(log_data_wo_burnin$B1)) {
      esses$B1 <- 400 
    }
    
  }
  
  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    next
    
  }else{
    inc(nr_converged_chains)
  }
  
  
  ### for the B1 index
  B1_index_hpd  = hdi(object = log_data_wo_burnin$'B1', credMass = 0.95)
  
  ### extract simulation value
  true_B1 = true_tree_stats[true_tree_stats$seed == seed, "true_b1"]

  ## check recovery
  recovered = is_parameter_in_hpd(hpd_lower = B1_index_hpd ["lower"], hpd_upper = B1_index_hpd ["upper"],
                                  true_parameter = true_B1)
  ## fill in table
  B1_index_inference[seed, ] = c(seed, B1_index_hpd , median(log_data_wo_burnin$'B1'), true_B1, recovered)
  
  
}

B1_index_inference = B1_index_inference[order(B1_index_inference$true_value), ]
B1_index_inference$orderedSeed = 1:100

#calculate coverage
sum(B1_index_inference$recovered)/nr_converged_chains

write.csv(x = B1_index_inference, "inference_logs/summary/B1_index_inference.csv",quote = F, row.names = F)






