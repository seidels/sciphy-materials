## Script name: plot_figure_3
##
## Purpose of script: Plot a figure that summarises how SciPhy compares with UPGMA
##
##

## set output directory for the plot
pic_dir = "plots/"

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
source("plot_trees_in_2d_with_likelihoods.R")

#extract likelihood values corresponding to the subsample of trees:

#### preprocess log to extract the correct likelihood values to plot ####
sciphy_trees <- ape::read.nexus(file = "inference_output/clockPerTarget_sampling_DataSet1_3000000.trees")
sample_nr_tree <- as.numeric(unlist(strsplit(names(sciphy_trees),"_"))[seq(2,2*length(sciphy_trees),by=2)])

## Load log data and extract likelihood values corresponding to the matching sample nrs
subsample_log <- read.table("inference_output/clockPerTarget_sampling_DataSet1_3000000.log",header = TRUE)
step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- subsample_log$Sample[2] - subsample_log$Sample[1]

##check if the logging frequency is the same
step_tree == step_log

#check that all tree sample nrs and likelihood sample nrs match
which(subsampled_log$Sample != sample_nr_tree)

#extract likelihood values
tree_likelihood <- subsampled_log$likelihood

#plot the 2D mapped distances:
CI_treeDist <- read.csv( file = "inference_output/mapping_df_PI.csv")

tree_df_rf <- CI_treeDist
plot1_2 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot1_2 <-  plot1_2 + theme(text=element_text(size = text_size),panel.border = element_blank() , axis.line = element_line() , axis.ticks = element_blank(), axis.text = element_blank() )
plot1_3 <-  plot1_3 + theme(text=element_text(size = text_size),panel.border = element_blank() , axis.line = element_line() , axis.ticks = element_blank(), axis.text = element_blank() )
text_size = 7
plot1_4  <-  plot1_4 + theme(text=element_text(size = text_size),panel.border = element_blank() , axis.line = element_line() , axis.ticks = element_blank(), axis.text = element_blank(),legend.position = c(0.5,0.5) )

legend <- get_legend(plot1_4) 
text_size = 10
plot1_4 <-  plot1_4 + theme(text=element_text(size = text_size),panel.border = element_blank() , axis.line = element_line() , axis.ticks = element_blank(), axis.text = element_blank(),legend.position = c(0.5,0.5) )


row_1 <- cowplot::plot_grid(plot1_2+ theme(plot.title = element_blank()),plot1_3 + theme(plot.title = element_blank()),plot1_4+ theme(plot.title = element_blank(),legend.position = "none"),legend ,nrow = 1,align = "h",label_size=text_size,rel_widths = c(1,1,1,0.5),rel_heights = c(1,1,1,0.5))
ggsave(paste0(pic_dir,"row_1.pdf"),
       row_1, width = 14.28, height = 7.14, units = "cm", dpi = 300)

# ----------------------------------------
## Plot the LTT plot: UPGMA against SciPhy sample 
## ----------------------------------------
##
##read in SciPhy trees
trees_constant <- ape::read.nexus(file = "inference_output/clockPerTarget_sampling_DataSet1_3000000.trees")

##read in upgma tree
upgma <- ape::read.tree(file = "inference_output/UPGMAtree_1000.txt")

##read in log_constant
log_constant <- read.table("../figure_3/inference_output/combined_clockPerTarget_sampling_DataSet1.log", header = T)

#get the median posterior tree height
median_posterior_height_constant <- median(log_constant[,"treeHeight.t.alignment"])

##read in the MCC_constant tree
MCC_constant <- ape::read.nexus(file = "inference_output/MCC_median_heights_clockPerTarget_sampling_DataSet1_3000000.tree")

##read in the MCC_constant tree
MCC_upgma_sciphy <- ape::read.nexus(file = "inference_output/combined_1000_UPGMA_medianPosteriorHeight_estimateBranchLengths_infer_rho_sampling_MCC.tree")



##read in scaled upgma
upgma_scaled_constant <- upgma
upgma_height <- tree_height_calc(upgma)
upgma_scaled_constant$edge.length <- upgma_scaled_constant$edge.length * (median_posterior_height_constant/upgma_height)

#write the scaled upgma
#ape::write.tree(upgma_scaled_constant, "UPGMA_scaled_clock_per_target_sampling.txt")

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

#add "off the shelf" upgma 
ltt_upgma <- ltt.plot.coords(upgma)
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_upgma,rep(length(trees_constant)+1,nrow(ltt_upgma))))

##add the scaled upgma (with the height scaled with the median of the "full inference" sciphy analysis)
ltt_upgma_scaled_constant <- ltt.plot.coords(upgma_scaled_constant)
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_upgma_scaled_constant,rep(length(trees_constant)+2,nrow(ltt_upgma_scaled_constant))))

#MCC main analysis, full inference
ltt_MCC_constant <- ltt.plot.coords(MCC_constant)   
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_MCC_constant,rep(length(trees_constant)+3,nrow(ltt_MCC_constant))))

#MCC UPGMA scaled branch lengths with SCIPHY
ltt_mcc_upgma_sciphy <- ltt.plot.coords(MCC_upgma_sciphy)   
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_mcc_upgma_sciphy,rep(length(trees_constant)+4,nrow(ltt_mcc_upgma_sciphy))))
all_ltt_constant <- data.frame(all_ltt_constant)

#making times start at zero (instead of a negative timescale)
all_ltt_constant$time <- all_ltt_constant$time + 25

#create another column to colour by type  
all_ltt_constant$type <- c(rep("SciPhy",length(which(all_ltt_constant$number<=(length(trees_constant))))),rep("UPGMA",length(which(all_ltt_constant$number==(length(trees_constant)+1)))),
                           rep("UPGMA root",length(which(all_ltt_constant$number==(length(trees_constant)+2)))),
                           rep("SciPhy MCC",length(which(all_ltt_constant$number==(length(trees_constant)+3)))),
                           rep("UPGMA + SciPhy",length(which(all_ltt_constant$number==(length(trees_constant)+4)))))

#define colors
cols <- c("SciPhy" = "#56996E", "SciPhy MCC" = "black", "UPGMA" = "red","UPGMA root" = "#E07E5E", "UPGMA + SciPhy" = "yellow" )

ltt_all <-  ggplot(all_ltt_constant,aes(x=time,y=N,group=number,colour=type,alpha = type)) + 
  geom_step() + 
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

## Add the UPGMA based es timates
typewriter_file <- "inference_output/combined_fixed_UPGMA_rho_sampling.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
HDInterval::hdi(bd_rates)
bd_rates$tree = "UPGMA"
bd_rates = melt(bd_rates, id.vars = "tree")
growth_rates = bd_rates

## Add UPGMA scaled based estimates
typewriter_file <- "inference_output/combined_1000_UPGMA_infer_rho_sampling_median_height.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA root"
bd_rates = melt(bd_rates, id.vars = "tree")
growth_rates = rbind(growth_rates, bd_rates)

## Add the main SciPhy analysis 
typewriter_file <- "../figure_3/inference_output/combined_clockPerTarget_sampling_DataSet1.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "SciPhy"
bd_rates = melt(bd_rates, id.vars = "tree")
growth_rates = rbind(growth_rates, bd_rates)

## Add the UPGMA with scaled branches with SciPhy
typewriter_file <- "inference_output/combined_1000_UPGMA_medianPosteriorHeight_estimateBranchLengths_infer_rho_sampling.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA +"
bd_rates = melt(bd_rates, id.vars = "tree")
growth_rates = rbind(growth_rates, bd_rates)


## Add the Prior
prior = data.frame(growthRate=rlnorm(n = 1000, meanlog = 0.1, sdlog = 1) -  rlnorm(1000, meanlog = -0.4, sdlog = 1))
prior$tree = "Prior"
prior =melt(prior, id.vars = "tree")
growth_rates = rbind(growth_rates, prior)

# Define order
growth_rates$tree = factor(growth_rates$tree, levels=c("Prior","SciPhy","Fixed sampling", "UPGMA", "UPGMA root","UPGMA +","SciPhy MCC"))
#define colors
cols <- c("Prior" = "white","Fixed sampling"="#B7CC62","SciPhy" = "#56996E", "UPGMA" = "red","UPGMA root" = "#E07E5E","SciPhy MCC" = "grey","UPGMA +" = "yellow","UPGMA root" = "pink","SciPhy MCC"="grey")

row_3 <- ggplot(growth_rates, aes(x=tree, y=value, fill = tree)) +
  theme_bw() +
  geom_violin(draw_quantiles = c(0.5))+  
  scale_fill_manual(values = cols) + 
  scale_y_continuous(breaks=c(0.2,0.6,1.0,1.4),limits=c(0.2, 1.6)) + 
  xlab("")+
  ylab(expression("Growth rate [" * d^-1 * "]")) + theme(legend.title=element_blank(),axis.text.x = element_text(color = "black", angle= 90, vjust = 0.7,hjust=1) )

ggsave(paste0(pic_dir,"row_3.pdf"),p_growth + theme(legend.position = "none") , width = 14.28, height = 7.14, units = "cm", dpi = 300)


####------------------------------------
#
# Combine and save all panels with labels
#
####-------------------------------------

full_figure <- cowplot::plot_grid(row_1,row_2,p_growth +theme(legend.position = "none"), ncol=1, labels = "AUTO", label_y=c(1,1.1,1.1),rel_heights = c(1,0.9,1.5)) 

ggsave(paste0(pic_dir,"figure_4.pdf"),full_figure, width = 14.28, height = 14.28, units = "cm", dpi = 300)




