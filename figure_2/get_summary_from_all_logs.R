## ---------------------------
##
## Script name: get summary data from all log files
##
## Purpose of script: Analyse well-calibrated simulations
##
## Author: Sophie Seidel
##
## Date Created: 2023-04-05
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##
## ---------------------------


library(tracerer)
library(HDInterval)
library(ggplot2)

## ---------------------------
log_dir = "inference_logs/"

simulation_dir = "simulation_parameters/"

parameters_of_interest = c("clockRate", paste0("insertRate.", 1:13))


inc = function(x){
  eval.parent(substitute(x <- x+1))
}

is_parameter_in_hpd = function(hpd_lower, hpd_upper, true_parameter){
  if(true_parameter >= hpd_lower && true_parameter <= hpd_upper)
    return(TRUE)
  else
    return(FALSE)
}

nr_converged_chains = 0

clock_rate_array_5_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

insert_rate_inference = data.frame(seed=rep(1:100, each = 13), hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,
                                   insertRate= rep(1:13,100))


for (seed in  1:100){

  print(seed)
  # get inference log
  log_file = paste0("infer_given_fixed_tree_13_inserts.", seed, ".log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.1)

  # get true parameter from simulation
  simulation_parameter_file = paste0(simulation_dir, "simParams_", seed, ".csv")
  simulation_parameters = read.csv(simulation_parameter_file)
  colnames(simulation_parameters) = c("parameter", "value")

  #check that chains has converged
  esses = calc_esses(log_data_wo_burnin, sample_interval = 1000)
  esses = esses[!colnames(esses) %in% c("birthRate", "deathRate", "prior")]
  esses = esses[! names(esses) %in% c("samplingProportion", "treeHeight.t.alignment", "treeLength.t.alignment")]
  esses

  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    next

  }else{
    inc(nr_converged_chains)
  }

  clock_rate_array_5_hpd = hdi(object = log_data_wo_burnin$clockRate, credMass = 0.95)

  true_clock_rate = simulation_parameters[14, 2]

  ## check for clock rate for alignments with array length 5
    recovered = is_parameter_in_hpd(hpd_lower = clock_rate_array_5_hpd["lower"], hpd_upper = clock_rate_array_5_hpd["upper"],
                      true_parameter = true_clock_rate)
  clock_rate_array_5_inference[seed, ] = c(seed, clock_rate_array_5_hpd, median(log_data_wo_burnin$clockRate), true_clock_rate, recovered)

  for (insert_rate_nr in 1:13){

    insert_rate = paste0("insertRates.", insert_rate_nr)
    insert_rate_hpd = hdi(object = log_data_wo_burnin[, insert_rate], credMass = 0.95)
    median_insert_rate = median(log_data_wo_burnin[, insert_rate])
    true_insert_rate = simulation_parameters[insert_rate_nr, "value"]
    recovered = is_parameter_in_hpd(hpd_lower = insert_rate_hpd["lower"], hpd_upper = insert_rate_hpd["upper"],
                                    true_parameter = true_insert_rate)


    row_index = which(insert_rate_inference$seed == seed & insert_rate_inference$insertRate == insert_rate_nr)
    insert_rate_inference[row_index, ] = c(seed, as.numeric(insert_rate_hpd), median_insert_rate, true_insert_rate,
                                           recovered, insert_rate_nr)
  }
}

### clock rate of alignment with length 5 targets
clock_rate_array_5_inference = clock_rate_array_5_inference[order(clock_rate_array_5_inference$true_value), ]
clock_rate_array_5_inference$orderedSeed = 1:100

#coverages per insert rate
coverages_per_insert <- c()
for(i in 1:13) {
  coverage <- sum(insert_rate_inference[which(insert_rate_inference$insertRate == i),"recovered"])/nr_converged_chains
  coverages_per_insert <- c(coverages_per_insert,coverage)
}
coverages_per_insert
#0.97 0.95 0.93 0.92 0.97 0.95 0.94 0.92 0.91 0.93 0.99 0.94 0.95

#coverage for the clock rate
sum(clock_rate_array_5_inference$recovered)/nr_converged_chains
#0.96

## true values on x axis and estimated values on y axis

g = ggplot(insert_rate_inference, aes(x=true_value, y=median)) +
  geom_point()+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median and posterior interval")+
  theme_classic()+
  xlim(0, 0.5) + ylim(0, 0.5)+
  ggtitle("Insertion probabilities")+
  theme(plot.title = element_text(hjust=0.5))
g
ggsave(plot = g, "plots/inference_insert_rate.pdf",
       width = 12, height = 12, units = "cm")

g = ggplot(clock_rate_array_5_inference, aes(x=true_value, y=median)) +
  geom_point()+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median and posterior interval")+
  theme_classic()+
  xlim(0, 0.5) + ylim(0, 0.5)+
  ggtitle("Editing rate")+
  theme(plot.title = element_text(hjust=0.5))
g
ggsave(plot = g, "plots/inference_clock_array_5.pdf",
       width = 12, height = 12, units = "cm")

