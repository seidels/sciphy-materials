## ---------------------------
##
##
## Purpose of script: plot mcc tree with nodes annotated with posterior support
##
## Date Created: 2023-10-13
##
## ---------------------------


## load up the packages we will need:  (uncomment as required)
library(ggtree)
library(ggplot2)
library(treeio)

# load tree metric libraries
library(treebalance)
library(treestats)
library(ape)

text_size = 10
output_dir = "plots/"
## ---------------------------


## define input files
tree_file = "inference_output/4-mGASv2-skyline-ou-40K.1000Kresampled-3-seeds.CCD0.defaultHeights.txt"

## load input
tree = read.beast(file = tree_file)
dat = tree@data

# Add origin branch to match experiment duration
origin_length = 11 - round(max(dat$height), digits = 3)

# plot tree
p = ggtree(tree, root.position = origin_length, size=0.15, color="darkgrey",ladderize = TRUE) + 
  theme_tree2() +
  geom_rootedge(rootedge = origin_length,size=0.15) +
  geom_nodepoint(aes(size=posterior)) + 
  scale_size(range = c(0.1, 4), name = "Posterior support") + 
  vexpand(.1, 1) +  vexpand(.01, -1) + 
  theme(text = element_text(size=text_size),legend.position=c(0.8,0.8)) + 
  scale_x_continuous(breaks = c(0, 4,6,7,7.5 ,11))+ xlab("Time [d]")

p

ggsave(paste0(output_dir,"fig_5_e_final.png"), p, width=12.8, height = 11, units = "cm", dpi = 300)
ggsave(paste0(output_dir,"fig_5_e_final.pdf"), p, width=12.8, height = 11, units = "cm", dpi = 300)

