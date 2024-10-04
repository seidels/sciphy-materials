## Script name: plot_figure_4_D 
##
## Purpose of script: Plot the gastruloid growth, under the 3 change time and 2 change time analysis with Merle et al. as input
##

setwd("/Users/azwaans/typewriter_analysis/results/gastruloid/analyses/")

## set working directory 
pic_dir = "~/typewriter_analysis/paper_figures/figure_5/"

# set figure settings
text_size=10
label_size = 12

## load up the packages we will need:  (uncomment as required)
library(tidyverse)
library(data.table)
library(ape)
library(TreeDist)
library(gdata)
library(reshape2)
library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)
library(pammtools)
library(reshape2)
library(rjson)
library(jsonlite)


get_median_and_hpd = function(growth){
  
  HPD <- HPDinterval(growth)
  
  name = paste0("growthRate.", 1:ncol(growth))
  median <- as.numeric(sapply( data.frame(growth),median))
  up_bd <- as.numeric(HPD[,"upper"])
  low_bd <- as.numeric(HPD[,"lower"])
  
  return(data.frame(name=name, median=median, hpd_up = up_bd, hpd_low=low_bd))
}


#load the raw data obtained from Merle et. al, 2024
merle_raw <- read_json(path="Extended_Data_Fig_2b_raw_data.json", simplifyVector = TRUE)
merle_raw <- merle_raw$data

#transform cell counts to growth rates. There are 23 replicates, with each 5 time points measurements
# if N0 is the number of cells at time t0 and N1 is the number of cells at time t1 (t1>t0) 
# the growth rate r under exponential growth is rate r = (1/(t1-t0))*log(N1/N0), time intervals for this dataset t1 - t0 = 1 day for all intervals
rate <- c()
for(i in 0:22) { 
  for(j in 2:5) {
    print(i*5 + j)
    print(merle_raw$Well_ID[i*5 + j])
    rate <- c(rate,log(merle_raw$N[i*5 + j]/merle_raw$N[i*5 + j-1]))
  }
}

#the start of the Merle et al experiment corresponds to day 6 in the Choi et al, experiment (1 day before CHIR pulse), assign those times:
time <- rep(c(6.0,7.0,8.0,9.0),23)
merle_time <- data.frame(T=time,growth=rate)

merle_mean <- aggregate(merle_time[, 2], list(merle_time$T), mean)
merle_sd <- aggregate(merle_time[, 2], list(merle_time$T), sd)
low <- merle_mean$x - merle_sd$x
high <- merle_mean$x + merle_sd$x

merle_low <- c(low,low[4]) 
merle_high <- c(high,high[4]) 
merle_mean <- c(merle_mean$x,merle_mean$x[4])

merle_et_al_full <- data_frame(t=c(6.0,7.0,8.0,9.0,10),mean=merle_mean,low=merle_low,high=merle_high, tree="Merle et al.,2024")

#recaculating everything ignoring intermediate points between 7 and 10 days. 

#transform cell counts to growth rates. There are 23 replicates, with each 5 time points measurements
# if N0 is the number of cells at time t0 and N1 is the number of cells at time t1 (t1>t0) 
# the growth rate r under exponential growth is rate r = (1/(t1-t0))*log(N1/N0), time intervals for this dataset t1 - t0 = 1 day for all intervals
rate <- c()
for(i in 0:22) { 
  # only use j=2 and j=5 to bin
  for(j in c(2,5)) {
    print(i*5 + j)
    print(merle_raw$Well_ID[i*5 + j])
    rate <- c(rate,log(merle_raw$N[i*5 + j]/merle_raw$N[i*5 + j-1]))
  }
}

#the start of the Merle et al experiment corresponds to day 6 in the Choi et al, experiment (1 day before CHIR pulse), assign those times:
time <- rep(c(6.0,7.0),23)
merle_time <- data.frame(T=time,growth=rate)

merle_mean <- aggregate(merle_time[, 2], list(merle_time$T), mean)
merle_sd <- aggregate(merle_time[, 2], list(merle_time$T), sd)
low <- merle_mean$x - merle_sd$x
high <- merle_mean$x + merle_sd$x

merle_low <- c(low,low[2]) 
merle_high <- c(high,high[2]) 
merle_mean <- c(merle_mean$x,merle_mean$x[2])

merle_et_al_bins <- data_frame(t=c(6.0,7.0,10),mean=merle_mean,low=merle_low,high=merle_high, tree="Merle et al.,2024")



# ---------------------------------------------------------
#  plot the dynamic growth rates inferred with and without the sampling proportion for the cell culture datasets
# ---------------------------------------------------------



timeline_changes <- c(0.0,4, 7.5, 11)

#load the log file with 2 change points in the estimated rates
log_file <- "combined_4-mGASv2-skyline-ou_1000000.log"
typewriter <- read.table(log_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:3)] - typewriter_mcmc[,paste0("deathRate.",1:3)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[4, ] =  growth_hpd[3, ]# add last timepoint, st for time - 11 the same value remains
growth_hpd$t = timeline_changes
growth_hpd$tree = "SciPhy"
growth_combined =   growth_hpd


#set colors and fills manually
cols <- c("SciPhy" = "red", "3 change times" = "pink","Merle et al.,2024"="black", "mean Merle et al."="lightblue","binned Merle et al."="blue")
cols_fill<-  c("SciPhy" = "red", "3 change times" = "pink","Merle et al.,2024"="black","mean Merle et al."="lightblue","binned Merle et al."="blue")
p_growth_OU_2 <- ggplot(growth_combined) +  
  geom_step(size=1,aes(x=t, y=median, col=tree))+ 
  geom_stepribbon(aes(x=t,ymin = hpd_low, ymax=hpd_up, fill=tree, col=tree), linetype="dotted", alpha = 0.1) +
  geom_step(data=merle_et_al_bins,aes(x=t,y=mean,col=tree),size=1,linetype="dashed")+
  geom_stepribbon(data=merle_et_al_bins,aes(x=t,ymin=merle_low, ymax=merle_high, fill=tree, col=tree), linetype="dotted", alpha = 0.1) +  
  theme_bw() + scale_color_manual(values = cols) + scale_fill_manual(values = cols_fill) + 
  theme(legend.title = element_blank(),text = element_text(size = text_size),panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line()) +
  ylab(expression("Growth rate [" * d^-1 * "]"))+ scale_x_continuous(breaks=c(0,4,6,7,10,11)) +
  xlab("Time [d]") 

p_growth_OU_2 <- p_growth_OU_2 + theme(legend.title = element_blank(),text = element_text(size = text_size),panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),legend.position = c(0.85,0.7)) 

#p_growth_OU <- p_growth_OU + theme(legend.position = "none")+  guides(color = guide_legend(ncol = 1)) ,ncol=2,nrow=1, rel_widths = c(1.0,1.0),axis = "l")
ggsave(paste0(pic_dir,"figure_5_D.pdf"),p_growth_OU_2, width = 14.28, height = 7.14, units = "cm", dpi = 300)


###########################
#Supplemental figure with the unbinned Merle et al. data and all time bins under SciPhy
###########################


#specify the timeline (breakpoints in for the skyline plot)
timeline_changes <- c(0.0,4, 7,8, 11)

#load the log file with 2 change points in the estimated rates

log_file <- "combined_3-mGASv2-skyline-ou_1000000.log"
typewriter <- read.table(log_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:4)] - typewriter_mcmc[,paste0("deathRate.",1:4)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[5, ] =  growth_hpd[4, ]# add last timepoint, st for time - 11 the same value remains
growth_hpd$t = timeline_changes
growth_hpd$tree = "3 change times"
growth_combined =  growth_hpd

timeline_changes <- c(0.0,4, 7.5, 11)

#load the log file with 2 change points in the estimated rates

log_file <- "combined_4-mGASv2-skyline-ou_1000000.log"
typewriter <- read.table(log_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:3)] - typewriter_mcmc[,paste0("deathRate.",1:3)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[4, ] =  growth_hpd[3, ]# add last timepoint, st for time - 11 the same value remains
growth_hpd$t = timeline_changes
growth_hpd$tree = "2 change times"
growth_combined = rbind(growth_combined,growth_hpd)


#set colors and fills manually
cols <- c("2 change times" = "red", "3 change times" = "orange","Merle et al.,2024"="black", "mean Merle et al."="lightblue","binned Merle et al."="blue")
cols_fill<-  c("2 change times" = "red", "3 change times" = "orange","Merle et al.,2024"="black","mean Merle et al."="lightblue","binned Merle et al."="blue")
p_growth_OU_2 <- ggplot(growth_combined) +  
  geom_step(size=1,aes(x=t, y=median, col=tree))+ 
  geom_stepribbon(aes(x=t,ymin = hpd_low, ymax=hpd_up, fill=tree, col=tree), linetype="dotted", alpha = 0.1) +
  geom_step(data=merle_et_al_full,aes(x=t,y=mean,col=tree),size=1)+
  geom_stepribbon(data=merle_et_al_full,aes(x=t,ymin=merle_low, ymax=merle_high, fill=tree, col=tree), linetype="dotted", alpha = 0.1) +  
  theme_bw() + scale_color_manual(values = cols) + scale_fill_manual(values = cols_fill) + 
  theme(legend.title = element_blank(),text = element_text(size = text_size),panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line()) +
  ylab(expression("Growth rate [" * d^-1 * "]"))+ scale_x_continuous(breaks=c(0,4,6,7,8,9,10,11)) +
  xlab("Time [d]") 

p_growth_OU_2 <- p_growth_OU_2 + theme(legend.title = element_blank(),text = element_text(size = text_size),panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),legend.position = c(0.85,0.7)) 

ggsave(paste0(pic_dir,"figure_5_supplement.pdf"),p_growth_OU_2, width = 14.28, height = 7.14, units = "cm", dpi = 300)


###
