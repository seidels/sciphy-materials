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
tree_file = "inference_output/4-mGASv2-skyline-ou.mcc.tree"

## load input
tree = read.beast(file = tree_file)
dat = tree@data

# Add origin branch to match experiment duration
origin_length = 11 - round(max(dat$height), digits = 3)

# plot tree
p = ggtree(tree, root.position = origin_length, size=0.05)  +
  geom_rootedge(rootedge = origin_length - 0.001, size=0.05) + 
  theme_tree2() + 
  theme(axis.text.x  = element_text(size = text_size),
        axis.title.x  = element_text(size = text_size),
        legend.position = c(0.1, 0.7) )+
  geom_nodepoint(aes(size=posterior))+
  scale_size(range = c(0.1, 2))+
  scale_x_continuous( breaks = c(0, 4, 7, 8, 11))+
  xlab("Time [d]")
p

ggsave(paste0(output_dir,"fig_5_e.png"), p, width = 14.28, height = 8, units = "cm", dpi = 300)
