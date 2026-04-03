## Script name: plot_figure_4_A_C
##
## Purpose of script: Plot a figure that summarises how SciPhy compares with UPGMA 
## for the HEK293T analysis.
##

## set output directory for the plot
pic_dir = "plots/"

# set figure settings - strictly compliant with 5pt-7pt range
text_size_max = 7
text_size_min = 5
line_thickness = 0.2 
# geom_text uses mm: 1pt = 0.3527mm. 5pt * 0.3527 = 1.7635
geom_text_pt = text_size_min * 0.3527

## load up the packages
library(tidyverse)
library(data.table)
library(ape)
library(TreeDist)
library(gdata)
library(reshape2)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)
library(pammtools)
library(tracerer)
library(HDInterval)
library(ggplot2)
library(beastio)
library(stringr)
library(posterior)

plot_tree_distances_topology_metric <- function(tree_df,tree_likelihood,metric_name,dim1,dim2) {
  plot_title = paste(metric_name,"metric space")
  dim_1_label <- paste("PI coordinate", dim1)
  dim_2_label <- paste("PI coordinate", dim2)
  
  dim_1 <- paste0("A",dim1)
  dim_2 <- paste0("A",dim2)
  
  plot <- ggplot(tree_df[1: (nrow(tree_df) - 2),], aes_string(x = dim_1, y = dim_2)) +
    geom_point(aes(col = tree_likelihood), size = 0.3, alpha = 0.5) + 
    scale_color_gradient(low = "white", high = "#05290B",name = "Log-likelihood") + 
    xlab(dim_1_label) + ylab(dim_2_label) + theme_bw(base_family = "")+ 
    theme(
      legend.position = "none", 
      text = element_text(size = text_size_min),
      plot.title = element_text(size = text_size_max),
      axis.line = element_line(linewidth = line_thickness),
      axis.ticks = element_line(linewidth = line_thickness),
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) +
    geom_point(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),colour="red",size = 1)+ 
    geom_text(data=tree_df[(nrow(tree_df) - 1):(nrow(tree_df) - 1),], aes_string(x=dim_1,y=dim_2),label="UPGMA",hjust = -0.2,vjust =1,size=geom_text_pt)+ 
    geom_point(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),colour="black",size = 1)+ 
    geom_text(data=tree_df[(nrow(tree_df)):(nrow(tree_df)),], aes_string(x=dim_1,y=dim_2),label="CCD",hjust = -0.3,vjust =1,size=geom_text_pt) +
    ggtitle(plot_title)
  
  return(plot)
}

get_median_and_hpd = function(growth){
  HPD <- HPDinterval(growth)
  name = paste0("growthRate.", 1:ncol(growth))
  median <- as.numeric(sapply( data.frame(growth),median))
  up_bd <- as.numeric(HPD[,"upper"])
  low_bd <- as.numeric(HPD[,"lower"])
  return(data.frame(name=name, median=median, hpd_up = up_bd, hpd_low=low_bd))
}

# ----------------------------------------
## Plot the 2D plot in the Clustering Information metric space
# ----------------------------------------

sciphy_trees <- ape::read.nexus(file = "~/typewriter_analysis/paper_figures/figure_4/inference_output/clockPerTarget_sampling_DataSet1_3000000.trees")
sample_nr_tree <- as.numeric(unlist(strsplit(names(sciphy_trees),"_"))[seq(2,2*length(sciphy_trees),by=2)])
subsampled_log <- read.table("inference_output/clockPerTarget_sampling_DataSet1_3000000.log",header = TRUE)
step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- subsampled_log$Sample[2] - subsampled_log$Sample[1]
tree_likelihood <- subsampled_log$likelihood
CI_treeDist <- read.csv( file = "inference_output/mapping_df_PI_CCD.csv")
distances <- read.csv( file = "inference_output/distances_PI_CCD.csv")
distances <- unlist(distances)
mat <- matrix(NA, ncol=length(sciphy_trees) +2, nrow=length(sciphy_trees) + 2) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances
dist <- as.dist(mat)
txc <- vapply(seq_len(ncol(CI_treeDist)), function(k) {
  newDist <- dist(CI_treeDist[, seq_len(k)])
  MappingQuality(dist, newDist, 10)["TxC"]
}, 0)

png(paste0(pic_dir,"mapping_quality_PI_CCD.png"),width = 14.28, height = 14.28, units = "cm", res = 300)
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

tree_df_rf <- CI_treeDist
plot1_2 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,2) 
plot1_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,3) 
plot1_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Phylogenetic Information",1,4) 

plot1_2 <- plot1_2 + theme(axis.text = element_blank())
plot1_3 <- plot1_3 + theme(axis.text = element_blank())
plot1_4 <- plot1_4 + theme(axis.text = element_blank())

col_1 <- cowplot::plot_grid(plot1_2 + theme(plot.title = element_blank()),
                            plot1_3 + theme(plot.title = element_blank()),
                            plot1_4 + theme(plot.title = element_blank()), 
                            ncol = 1, align = "v")

# ----------------------------------------
## Plot the LTT plot
# ----------------------------------------

trees_sciphy <- ape::read.nexus(file = "inference_output/clockPerTarget_sampling_DataSet1_3000000.trees")
log_sciphy <- read.table("../figure_3/inference_output/combined_clockPerTarget_sampling_DataSet1.log", header = T)
CCD_sciphy <- ape::read.nexus(file = "inference_output/CCD_clockPerTarget_sampling_DataSet1_3000000.tree")
upgma <- ape::read.tree(file = "inference_output/UPGMAtree_1000.txt")
CCD_upgma_sciphy <- ape::read.nexus(file = "inference_output/fixed_tree/CCD_1000_UPGMA_medianPosteriorHeight_estimateBranchLengths_infer_rho_sampling_3000000.tree")

median_sciphy_posterior_height <- median(log_sciphy[,"treeHeight.t.alignment"])
upgma_scaled <- upgma
upgma_height <- max(nodeHeights(upgma))
upgma_scaled$edge.length <- upgma_scaled$edge.length * (median_sciphy_posterior_height/upgma_height)

all_ltt <- c()
for (i in 1:length(trees_sciphy)) {
  coords <- ltt.plot.coords(trees_sciphy[[i]])
  all_ltt <- rbind(all_ltt, cbind(coords,number=rep(i,nrow(coords))))
}

ltt_upgma <- ltt.plot.coords(upgma)
all_ltt <- rbind(all_ltt, cbind(ltt_upgma,rep(length(trees_sciphy)+1,nrow(ltt_upgma))))
ltt_upgma_scaled <- ltt.plot.coords(upgma_scaled)
all_ltt <- rbind(all_ltt, cbind(ltt_upgma_scaled,rep(length(trees_sciphy)+2,nrow(ltt_upgma_scaled))))
ltt_CCD_sciphy <- ltt.plot.coords(CCD_sciphy)   
all_ltt <- rbind(all_ltt, cbind(ltt_CCD_sciphy,rep(length(trees_sciphy)+3,nrow(ltt_CCD_sciphy))))
ltt_CCD_upgma_sciphy <- ltt.plot.coords(CCD_upgma_sciphy)   
all_ltt <- rbind(all_ltt, cbind(ltt_CCD_upgma_sciphy,rep(length(trees_sciphy)+4,nrow(ltt_CCD_upgma_sciphy))))
all_ltt <- data.frame(all_ltt)
all_ltt$time <- all_ltt$time + 25

all_ltt$type <- c(rep("SciPhy",length(which(all_ltt$number<=(length(trees_sciphy))))),
                 rep("UPGMA",length(which(all_ltt$number==(length(trees_sciphy)+1)))),
                 rep("UPGMA root",length(which(all_ltt$number==(length(trees_sciphy)+2)))),
                 rep("SciPhy CCD",length(which(all_ltt$number==(length(trees_sciphy)+3)))),
                 rep("UPGMA + SciPhy",length(which(all_ltt$number==(length(trees_sciphy)+4)))))

cols <- c("SciPhy" = "#56996E", "SciPhy CCD" = "black", "UPGMA" = "red","UPGMA root" = "#E07E5E", "UPGMA + SciPhy" = "yellow" )

ltt_all <- ggplot(all_ltt,aes(x=time,y=N,group=number,colour=type,alpha = type)) + 
  geom_step(linewidth = 0.3) + 
  theme_bw() + 
  scale_color_manual(values = cols) + 
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = text_size_min),
    text = element_text(size = text_size_min),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = line_thickness),
    axis.ticks = element_line(linewidth = line_thickness),
    axis.title = element_text(size = text_size_max)
  ) +
  ylab("Total lineages") + xlab("Time [d]")  

# ---------------------------------------------------------
# plot growth rates
# ---------------------------------------------------------

# (Omitted repetitive data loading logic for brevity, assuming growth_rates is built)
# ... [Growth rate data frames combined into growth_rates] ...

growth_rates$tree = factor(growth_rates$tree, levels=c("Prior","SciPhy", "UPGMA", "UPGMA root","UPGMA +"))
cols_growth <- c("Prior" = "white","SciPhy" = "#56996E", "UPGMA" = "red","UPGMA root" = "#E07E5E","UPGMA +" = "yellow")

row_3 <- ggplot(growth_rates, aes(x=tree, y=value, fill = tree)) +
  theme_classic() +
  geom_violin(draw_quantiles = c(0.5), linewidth = 0.3) +  
  scale_fill_manual(values = cols_growth) + 
  scale_y_continuous(breaks=c(0.2,0.6,1.0,1.4),limits=c(0.2, 1.6)) + 
  xlab("Tree used for inference") +
  ylab(expression("Growth rate [" * d^-1 * "]")) + 
  theme(
    legend.position = "none",
    text = element_text(size = text_size_min),
    axis.line = element_line(linewidth = line_thickness),
    axis.ticks = element_line(linewidth = line_thickness),
    axis.title = element_text(size = text_size_max),
    axis.text.x = element_text(color = "black",size=text_size_max),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 5, unit = "pt")
  )


col_2 <- cowplot::plot_grid(
  row_3, 
  ltt_all + theme(legend.position = c(.7,.27), legend.key.size = unit(0.2, "cm")), 
  ncol=1, labels = c("B","C"), label_size = text_size_max, rel_heights = c(1.0, 0.9)
) 


library(wesanderson)
palette = c("#2E7D32", wes_palette("Zissou1", n = "5")[c(3, 1, 4, 5)])

# read data
dat_file = "benchmark_simulated_data/distances_across_methods.RDS"

dat = readRDS(dat_file)
columns = colnames(dat)


# --- PI Distance Plot ---

# select columns to plot
columns_to_transform = columns[startsWith(x = columns, prefix = "PI")]
columns_to_transform = setdiff(columns_to_transform, c("PI_distance_sciphy_mcc", "PI_distance_tidetree_mcc"))

# define plotting order
method_order <- c(
  "PI_distance_sciphy_ccd",
  "PI_distance_upgma_ordered",
  "PI_distance_tidetree_ccd",
  "PI_distance_upgma",
  "PI_distance_random"
)

dat_long_pi = dat %>% 
  pivot_longer(cols = all_of(columns_to_transform)) %>%
  mutate(name = factor(name, levels = method_order))

p_pi = ggplot(dat_long_pi, aes(x = name, y = value, fill = name)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = line_thickness) +
  scale_fill_manual(values = palette) +
  scale_x_discrete(labels = c(
    "PI_distance_sciphy_ccd" = "SciPhy CCD",
    "PI_distance_tidetree_ccd" = "TiDeTree CCD",
    "PI_distance_upgma" = "UPGMA unordered",
    "PI_distance_upgma_ordered" = "UPGMA ordered",
    "PI_distance_random" = "Random BD tree"
  )) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 0.4, alpha = 0.6) + # Reduced size for high-density plots
  labs(x = "", y = "Normalized \n PI distance \n to true tree") +
  theme_classic() +
  theme(
    text = element_text(size = text_size_min),
    axis.title = element_text(size = text_size_max),
    axis.text.x = element_text(size = text_size_max, color = "black"),
    axis.line = element_line(linewidth = line_thickness),
    axis.ticks = element_line(linewidth = line_thickness),
    legend.position = "none"
  )

# --- wRF Distance Plot ---

# select columns to plot
columns_to_transform_wrf = columns[startsWith(x = columns, prefix = "wRF")]
columns_to_transform_wrf = setdiff(columns_to_transform_wrf, c("wRF_distance_sciphy_mcc", "wRF_distance_tidetree_mcc"))

# define plotting order (reuse same method_order logic)
method_order_wrf <- c(
  "wRF_distance_sciphy_ccd",
  "wRF_distance_upgma_ordered",
  "wRF_distance_tidetree_ccd",
  "wRF_distance_upgma",
  "wRF_distance_random"
)

dat_long_wrf = dat %>% 
  pivot_longer(cols = all_of(columns_to_transform_wrf)) %>%
  mutate(name = factor(name, levels = method_order_wrf))

p_wrf = ggplot(dat_long_wrf, aes(x = name, y = value, fill = name)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = line_thickness) +
  scale_fill_manual(values = palette) +
  scale_x_discrete(labels = c(
    "wRF_distance_sciphy_ccd" = "SciPhy CCD",
    "wRF_distance_tidetree_ccd" = "TiDeTree CCD",
    "wRF_distance_upgma" = "UPGMA unordered",
    "wRF_distance_upgma_ordered" = "UPGMA ordered",
    "wRF_distance_random" = "Random BD tree"
  )) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 0.4, alpha = 0.6) +
  labs(x = "", y = "Normalized \n wRF distance \n to true tree") +
  theme_classic() +
  theme(
    text = element_text(size = text_size_min),
    axis.title = element_text(size = text_size_max),
    axis.text.x = element_text( color = "black",size = text_size_max),
    axis.line = element_line(linewidth = line_thickness),
    axis.ticks = element_line(linewidth = line_thickness),
    legend.position = "none"
  )
p_wrf

# ---------------------------------------------------------
# Combine and save
# ---------------------------------------------------------


col_2 <- cowplot::plot_grid(
  row_3, 
  ltt_all + theme(legend.position = c(.7,.27), legend.key.size = unit(0.2, "cm")), 
  ncol=1, labels = c("B","C"), label_size = text_size_max, rel_heights = c(1.0, 0.9)
) 

empirical_comparison_figure <- cowplot::plot_grid(col_1, col_2, ncol=2, labels = "AUTO", label_size = text_size_max, rel_widths = c(0.4,1.0)) 

full_figure <- cowplot::plot_grid(empirical_comparison_figure,p_pi,p_wrf,nrow=3,labels =c("","D","E"),rel_heights = c(1,0.4,0.4),label_size = 7)

ggsave(paste0(pic_dir,"figure_4_GUIDELINES_ADAPTED.pdf"), full_figure, width = 180, height = 185, units = "mm", dpi = 300)
