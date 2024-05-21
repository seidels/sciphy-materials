## Script name: plot_figure_3
##
## Purpose of script: Plot a figure that summarises how SciPhy compares with UPGMA
##
##

## set working directory 
setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/TreeDist/")

## set output directory for the plot
pic_dir = "~/typewriter_analysis/paper_figures/figure_3/"

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


get_median_and_hpd = function(growth){
  
  HPD <- HPDinterval(growth)
  
  name = paste0("growthRate.", 1:ncol(growth))
  median <- as.numeric(sapply( data.frame(growth),median))
  up_bd <- as.numeric(HPD[,"upper"])
  low_bd <- as.numeric(HPD[,"lower"])
  
  return(data.frame(name=name, median=median, hpd_up = up_bd, hpd_low=low_bd))
}

# ----------------------------------------
## Plot the 2D plot in the Clustering Information metric stpace
## ----------------------------------------
##
##sourcing helper functions
source("~/typewriter_analysis/src/cell_culture/upgma/plot_trees_in_2d_with_likelihoods.R")

#extract likelihood values corresponding to the subsample of trees:

#### preprocess log to extract the correct likelihood values to plot ####
sciphy_trees <- ape::read.nexus(file = "~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/thinned4000000.trees")

sample_nr_tree <- as.numeric(unlist(strsplit(names(sciphy_trees),"_"))[seq(2,2*length(sciphy_trees),by=2)])

## Load log data and extract likelihood values corresponding to the matching sample nrs
log <- read.table("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/combined.log", header = T)

step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- log$Sample[2] - log$Sample[1]

#resample the log file at the same frequency
subsampled_log <- log[seq(1,length(log$Sample),by=step_tree/step_log),]

#resample the same number
subsampled_log <- subsampled_log[1:min(length(sciphy_trees),length(subsampled_log$Sample)),]

#check that all tree sample nrs and likelihood sample nrs match
which(subsampled_log$Sample != sample_nr_tree)

#extract likelihood values
tree_likelihood <- subsampled_log$likelihood

#plot the 2D mapped distances:

CI_treeDist <- read.csv( file = "mapping_df_CI.csv")

tree_df_rf <- CI_treeDist
plot1_2 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot1_2 <-  plot1_2 + theme(text=element_text(size = text_size),panel.border = element_blank() , axis.line = element_line() , axis.ticks = element_blank(), axis.text = element_blank() )
plot1_3 <-  plot1_3 + theme(text=element_text(size = text_size),panel.border = element_blank() , axis.line = element_line() , axis.ticks = element_blank(), axis.text = element_blank() )
plot1_4 <-  plot1_4 + theme(text=element_text(size = text_size),panel.border = element_blank() , axis.line = element_line() , axis.ticks = element_blank(), axis.text = element_blank(),legend.position = c(0.1,0.1) )

legend <- get_legend(plot1_4, position= "right") 

row_1 <- cowplot::plot_grid(plot1_2+ theme(plot.title = element_blank()),plot1_3 + theme(plot.title = element_blank()),plot1_4+ theme(plot.title = element_blank(),legend.position = "none"),legend,nrow = 1,align = "h",label_size=text_size,rel_widths = c(1,1,1,1))
ggsave(paste0(pic_dir,"row_1.pdf"),
       row_1, width = 14.28, height = 7.14, units = "cm", dpi = 300)

# ----------------------------------------
## Plot the LTT plot: UPGMA against SciPhy sample 
## ----------------------------------------
##
setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/")

##read in SciPhy trees
trees_constant <- ape::read.nexus(file = "thinned4000000.trees")

##read in upgma tree
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")

##read in log_constant
log_constant <- read.table("combined.log", header = T)

##read in the MCC_constant tree
MCC_constant <- ape::read.nexus(file = "MCC_on_thinned4000000_check.tree")

#get the median posterior tree height
median_posterior_height_constant <- median(log_constant[,"treeHeight.t.alignment"])

##read in scaled upgma
upgma_scaled_constant <- upgma
upgma_height <- tree_height_calc(upgma)
upgma_scaled_constant$edge.length <- upgma_scaled_constant$edge.length * (median_posterior_height_constant/upgma_height)

#extract ltt coordinates from all trees_constant
all_ltt_constant <- c()
#record at what time point each tree reaches 500 cells and 1000 cells
all_500_times <- c()
all_1000_times <- c()
for (i in 1:length(trees_constant)) {
  coords <- ltt.plot.coords(trees_constant[[i]])
  all_500_times <- c(all_500_times,coords[500,'time'])
  all_1000_times <- c(all_1000_times,coords[1000,'time'])
  all_ltt_constant <- rbind(all_ltt_constant, cbind(coords,number=rep(i,nrow(coords))))
}

#add upgma 
ltt_upgma <- ltt.plot.coords(upgma)
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_upgma,rep(length(trees_constant)+1,nrow(ltt_upgma))))

##add the scaled upgma
ltt_upgma_scaled_constant <- ltt.plot.coords(upgma_scaled_constant)
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_upgma_scaled_constant,rep(length(trees_constant)+2,nrow(ltt_upgma_scaled_constant))))

#MCC_constant
ltt_MCC_constant <- ltt.plot.coords(MCC_constant)   
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_MCC_constant,rep(length(trees_constant)+3,nrow(ltt_MCC_constant))))

all_ltt_constant <- data.frame(all_ltt_constant)

#making times start at zero (instead of a negative timescale)
all_ltt_constant$time <- all_ltt_constant$time + 25

#create another column to colour by type  
all_ltt_constant$type <- c(rep("SciPhy constant",length(which(all_ltt_constant$number<=(length(trees_constant))))),rep("UPGMA origin",length(which(all_ltt_constant$number==(length(trees_constant)+1)))),rep("UPGMA root (constant)",length(which(all_ltt_constant$number==(length(trees_constant)+2)))),rep("SciPhy MCC",length(which(all_ltt_constant$number==(length(trees_constant)+3)))))
cols <- c("SciPhy constant" = "#E1E1F7", "SciPhy MCC" = "black", "UPGMA origin" = "red","UPGMA root (constant)" = "#E07E5E")


## set working directory to where trees and log_constant files are

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target_OU/")    

##read in SciPhy trees
trees <- ape::read.nexus(file = "thinned_5000000.trees")

##read in upgma tree
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")

##read in log
log <- read.table("combined_burnin10.log", header = T)

##read in the mcc tree
MCC <- ape::read.nexus(file = "MCC_median_thinned_5000000.tree")

#get the median posterior tree height
median_posterior_height <- median(log[,"treeHeight.t.alignment"])

##read in scaled upgma
upgma_scaled <- upgma
upgma_height <- tree_height_calc(upgma)
upgma_scaled$edge.length <- upgma_scaled$edge.length * (median_posterior_height/upgma_height)

#extract ltt coordinates from all trees
all_ltt <- c()
#record at what time point each tree reaches 500 cells and 1000 cells
all_500_times <- c()
all_1000_times <- c()
for (i in 1:length(trees)) {
  coords <- ltt.plot.coords(trees[[i]])
  all_500_times <- c(all_500_times,coords[500,'time'])
  all_1000_times <- c(all_1000_times,coords[1000,'time'])
  all_ltt <- rbind(all_ltt, cbind(coords,number=rep(i,nrow(coords))))
}

#add upgma 
ltt_upgma <- ltt.plot.coords(upgma)
all_ltt <- rbind(all_ltt, cbind(ltt_upgma,rep(length(trees)+1,nrow(ltt_upgma))))

##add the scaled upgma
ltt_upgma_scaled <- ltt.plot.coords(upgma_scaled)
all_ltt <- rbind(all_ltt, cbind(ltt_upgma_scaled,rep(length(trees)+2,nrow(ltt_upgma_scaled))))

#MCC
ltt_mcc <- ltt.plot.coords(MCC)   
all_ltt <- rbind(all_ltt, cbind(ltt_mcc,rep(length(trees)+3,nrow(ltt_mcc))))

all_ltt <- data.frame(all_ltt)

#making times start at zero (instead of a negative timescale)
all_ltt$time <- all_ltt$time + 25


all_ltt$type <- c(rep("SciPhy 13-epochs",length(which(all_ltt$number<=(length(trees))))),rep("UPGMA origin",length(which(all_ltt$number==(length(trees)+1)))),rep("UPGMA root (13-epochs)",length(which(all_ltt$number==(length(trees)+2)))),rep("SciPhy MCC",length(which(all_ltt$number==(length(trees)+3)))))
cols <- c("SciPhy constant" = "#E1E1F7", "UPGMA origin" = "red","UPGMA root (constant)" = "#E07E5E","SciPhy MCC" = "black", "SciPhy 13-epochs" = "#56996E", "UPGMA" = "red","UPGMA root (13-epochs)" = "#E07E5E")
alphas <- c("SciPhy constant" = 0.3, "UPGMA origin" = 1.0,"UPGMA root (constant)" = 1.0,"SciPhy MCC" = 1.0,"SciPhy 13-epochs" = 0.3, "UPGMA" = 1.0,"UPGMA root (13-epochs)" = 1.0)

ltt_all <-  ggplot(all_ltt_constant,aes(x=time,y=N,group=number,colour=type,alpha = type)) + 
  geom_step() + 
  geom_step(data=all_ltt) +
  theme_bw() + 
  scale_color_manual(values = cols) + 
  theme(legend.title = element_blank(),text = element_text(size = text_size),panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line()) +
  ylab("Total lineages") + 
   xlab("Time [d]") + scale_alpha_manual(values=alphas) 

legend_1 <- get_legend(ltt_all) 
row_2 <- cowplot::plot_grid(ltt_all + theme(legend.position = "none"),legend_1,ncol=2,nrow=1, rel_widths = c(1.0,0.4))
ggsave(paste0(pic_dir,"row_2.pdf"), row_2, width = 14.28, height = 7.14, units = "cm", dpi = 300)



# ---------------------------------------------------------
#  plot the constant growth rates inferred on fixed trees
# ---------------------------------------------------------

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/fixed_tree_analyses/fixed_trees")

## Add the UPGMA based estimates
typewriter_file <- "1000_UPGMA.1707411001885.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA origin"

growth_rates = melt(bd_rates, id.vars = "tree")

## Add UPGMA scaled based estimates
typewriter_file <- "1000_UPGMA_medianPosteriorHeight.1707410740300.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA root"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the MCC tree based estimates
typewriter_file <- "1000_MCC.1707412619101.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "SciPhy MCC"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the main analysis estimates
typewriter_file <- "../../inference_results/clock_per_target/1000_cells/combined.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "SciPhy posterior"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the Prior
prior = data.frame(growthRate=rlnorm(n = 1000, meanlog = -0.6, sdlog = 1))
prior$tree = "Prior"
prior =melt(prior, id.vars = "tree")

growth_rates = rbind(growth_rates, prior)

# Define order
growth_rates$tree = factor(growth_rates$tree, levels=c("Prior", "SciPhy posterior", "SciPhy MCC", "UPGMA origin", "UPGMA root"))
#define colors
cols <- c("Prior" = "white","SciPhy posterior" = "#E1E1F7", "UPGMA origin" = "red","UPGMA root" = "#E07E5E","SciPhy MCC" = "grey")
p_growth <- ggplot(growth_rates, aes(x=tree, y=value, fill = tree)) +
  theme_bw() +
  geom_violin(draw_quantiles = c(0.5))+  
  scale_fill_manual(values = cols) + 
  ylim(0.3, 0.5)+
  xlab("")+
  ylab(expression("Constant growth rate [" * d^-1 * "]")) + theme(legend.title=element_blank())
 

# ---------------------------------------------------------
#  plot the dynamic growth rates inferred on fixed trees
# ---------------------------------------------------------

#specify the timeline (breakpoints in for the skyline plot)
timeline_changes <- seq(0,25,by=25/13)

#load the log file with 3 change points in the estimated rates
log_file <- "1000_MCC.combined.log"
typewriter <- read.table(log_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] =  growth_hpd[13, ]# add last timepoint, st for time - 11 the same value remains
growth_hpd$t = timeline_changes
growth_hpd$tree = "SciPhy MCC"

growth_combined = growth_hpd

#load the log file with 3 change points in the estimated rates
log_file <- "1000_UPGMA_medianPosteriorHeight.combined.log"
typewriter <- read.table(log_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] =  growth_hpd[13, ]# add last timepoint, st for time - 11 the same value remains
growth_hpd$t = timeline_changes
growth_hpd$tree = "UPGMA root"

growth_combined = rbind(growth_combined, growth_hpd)


#load the log file with 2 change points in the estimated rates
log_file <- "1000_UPGMA.combined.log"
typewriter <- read.table(log_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] =  growth_hpd[13, ]# add last timepoint, st for time - 11 the same value remains
growth_hpd$t = timeline_changes
growth_hpd$tree = "UPGMA origin"
growth_combined = rbind(growth_combined, growth_hpd)

#set colors and fills manually
cols <- c("Prior" = "white","SciPhy posterior" = "#E1E1F7", "UPGMA origin" = "red","UPGMA root" = "#E07E5E","SciPhy MCC" = "grey")
cols_fill<- c("Prior" = "white","SciPhy posterior" = "lightgreen", "UPGMA origin" = "red","UPGMA root" = "brown","SciPhy MCC" = "grey")

p_growth_OU = ggplot(growth_combined, aes(x=t, y=median, col=tree)) +
  geom_step(size=1)+
  geom_stepribbon(aes(ymin = hpd_low, ymax=hpd_up, fill=tree), linetype="dotted", alpha = 0.1)+
  theme_bw() + scale_color_manual(values = cols) + scale_color_manual(values = cols_fill) + 
  theme(legend.title=element_blank(), legend.position = c(0.27, 0.858)) +
  ylab(expression("Dynamic growth rate [" * d^-1 * "]"))+
  xlab("Time [d]")

legend <- get_legend(p_growth) 
p_growth_OU <- p_growth_OU + theme(text=element_text(size = text_size),panel.border = element_blank() , axis.line = element_line() ,panel.grid = element_blank())

row_3 <- cowplot::plot_grid(p_growth+ theme(legend.position = "none") ,p_growth_OU + theme(legend.position = "none"),legend,ncol=3,nrow=1, rel_widths = c(1.0,1.0,0.6),axis = "l")

ggsave(paste0(pic_dir,"row_3.pdf"),full_figure, width = 14.28, height = 7.14, units = "cm", dpi = 300)

####------------------------------------
#
# Combine and save all panels with labels
#
####-------------------------------------

full_figure <- cowplot::plot_grid(row_1,row_2,row_3, ncol=1, labels = "AUTO") 

ggsave(paste0(pic_dir,"figure_3.pdf"),full_figure, width = 14.28, height = 14.28, units = "cm", dpi = 300)




