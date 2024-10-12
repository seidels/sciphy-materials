## ---------------------------
##
## Purpose of script: contrast growth rate estimated with and 
## without sequence data
##
## Author: Sophie Seidel
##
## Date Created: 2023-10-15
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##

# set figure settings
text_size=10
figure_path = "plots/supp_fig_11.pdf"

# load inference logs
log_seq_data <- "../figure_3/inference_output/1-combined.log"
dat_seq_data <- read.table(log_seq_data, header = T)
dat_seq_data$growth_rate = dat_seq_data$birthRate - dat_seq_data$deathRate
dat_seq_data$analysis = "with seq data"

log_no_seq_data <- "supp_fig_11_data/typewriter_clockPerSite_13Sites_1000Cells_DataSet1_samplingFromPrior.1.log"
dat_no_seq_data <- read.table(log_no_seq_data, header = T)
dat_no_seq_data$growth_rate = dat_no_seq_data$birthRate - dat_no_seq_data$deathRate
dat_no_seq_data$analysis = "w/o seq data"

dat_combined = rbind(dat_seq_data[, c("growth_rate", "analysis")], dat_no_seq_data[, c("growth_rate", "analysis")])

growth_plots = 
  ggplot(dat=dat_combined, aes(x=analysis, y=growth_rate))+
  geom_violin()+
  ylim(0.4, 0.87)+
  ylab("Posterior growth rates per day")+
  theme_bw()

ggsave(filename = figure_path, plot = growth_plots)
