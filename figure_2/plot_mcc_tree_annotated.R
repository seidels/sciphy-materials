## ---------------------------
##
## Script name: plot mcc_tree_annotated
##
## Purpose of script: plot tree facing the typewriter sequences, with barcode annotations and colors matching insertion probability plots
##
## Date Created: 2023-10-13
##
## ---------------------------
##
## Notes: based on plot_mcc_tree.R
##
## ---------------------------


## load up the packages we will need:  (uncomment as required)
library(ggtree)
library(ggplot2)
library(treeio)

text_size = 10
output_dir = "~/Projects/sciphy-materials/figure_2/"
## ---------------------------

swap_integer_for_edit = function(integer, insert_to_integer_map){

  insert = insert_to_integer_map[insert_to_integer_map$integer == integer, "insert"]
  return(insert)
}

## ---------------------------

## define input files

tree_file = "~/Projects/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/MCC_on_thinned4000000trees.tree"
alignment = "~/Projects/typewriter_analysis/results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/alignment_seed1.txt"
cell_ids_file = "~/Projects/typewriter_analysis/results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/cell_ids_seed1.txt"
edit_file = "~/Projects/typewriter_analysis/results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/edit_table_sample_1.csv"
edit_to_integer_map = "~/Projects/typewriter_analysis/data/cell_culture/insert_to_integer_map.csv"


# Hack from github due to package error: could not find function "offspring.tbl_tree_item"
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

## load input
tree = read.beast(file = tree_file)
dat = tree@data

# show overview figure
# p = ggtree(tree)  +
#   geom_hilight(node=1004, fill="steelblue", alpha=0.5)+
#   geom_rootedge(rootedge = 0.35223) + theme_tree2()
# p
# 
# p = ggtree(tree, root.position = 0.35223)  +
#   geom_hilight(node=1962, fill="steelblue", alpha=0.5)+
#   geom_rootedge(rootedge = 0.3522) + theme_tree2() + geom_nodelab(aes(x=branch, label=node))
# p

p = ggtree(tree, root.position = 0.35223, size=0.05)  +
  geom_hilight(node=1478, fill="steelblue", alpha=0.5)+
  geom_rootedge(rootedge = 0.3522, size=0.05) + theme_tree2() + 
  theme(axis.text.x  = element_text(size = text_size-5),
        axis.title.x  = element_text(size = text_size))+
  xlab("Time [d]")
# + geom_nodelab(aes(x=branch, label=node))
p

ggsave(paste0(output_dir,"MCC_tree_marked_clade.png"), p, width = 4.76, height = 4.76, units = "cm", dpi = 300)
svg(filename = paste0(output_dir,"MCC_stree_marked_clade.svg"), width = 14.28, height = 15)
p
dev.off()

# zoom in on clade and show alignment and time uncertainty
sub =  tree_subset(tree = tree, node = 1962,  group_node = T, root_edge = T, levels_back = F)
sub =  tree_subset(tree = tree, node = 1478,  group_node = T, root_edge = T, levels_back = F)
subdat = sub@data
subtree = get.tree(sub)

# sub tree with negative branch lengths
sub_tree_plot = ggtree(sub, root.position = 3.8) + theme_tree2() + geom_rootedge(rootedge = 3.8) +  
  geom_range('height_0.95_HPD', color='grey', size=3, alpha=.4)+
  geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.1, size=3) 

tips_in_subtree = subtree$tip.label

# get alignment
insert_to_integer_map = read.csv(edit_to_integer_map)

cell_ids = read.csv(cell_ids_file, header = F)
cell_ids$tip_label = 0:999 # keep numbering of tips as used in alignment

# get cells that are present in subtree
cell_ids = cell_ids[which(cell_ids$tip_label %in% tips_in_subtree), ]
cell_ids = cell_ids[order(cell_ids$V1), ] # order, s.t. edits can be combined relying on ordered cell names


edits = read.csv(edit_file)

## create targetbc matrices
ctr = 1
targetBCs = unique(edits$TargetBC)
targetbc_edits_list =  vector(mode = "list", length = length(targetBCs))

for (targetBC in targetBCs){

  print(targetBC)
  edits_tbc = edits[ edits$TargetBC == targetBC, ]
  edits_tbc = edits_tbc[which(edits_tbc$Cell %in% cell_ids$V1), ]

  # add cell ids ;
  ## order cell barcodes alphabetically and the reuse the tip labels from cell ids
  edits_tbc = edits_tbc[order(edits_tbc$Cell), ]
  rownames(edits_tbc) = cell_ids$tip_label
  edits_tbc = edits_tbc[ , 4:8]

  # convert edits from integers to trinucleotides
  edits_tbc_trinucl = data.frame(Site1 = sapply(edits_tbc$Site1, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
                             Site2 = sapply(edits_tbc$Site2, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
                             Site3 = sapply(edits_tbc$Site3, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
                             Site4 = sapply(edits_tbc$Site4, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
                             Site5 = sapply(edits_tbc$Site5, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)})
  )
  rownames(edits_tbc_trinucl) = rownames(edits_tbc)

  targetbc_edits_list[[ctr]] = edits_tbc_trinucl
  ctr = ctr + 1
}

## plot tree with edit matrix
#basic_tree = ggtree(sub,layout="ellipse")

sub_tree_plot = ggtree(sub, root.position = 3.8) + theme_tree2() + geom_rootedge(rootedge = 3.8) +  
  geom_range('height_0.95_HPD', color='grey', size=3, alpha=.4)+
  geom_nodepoint(aes(size=posterior))+ scale_size(range = c(0.1, 2))
  

sub_tree_plot

for (i in 1:length(targetBCs)){
  print(i)
  labels <- c(targetBCs[i],rep("",ncol(targetbc_edits_list[[i]])-1))
  annotated_basic_tree = gheatmap(sub_tree_plot, targetbc_edits_list[[i]],
                                  hjust=0,colnames_offset_y=-10,colnames_offset_x=-10,
                                  colnames_position = "bottom",width = 0.05, offset = 1.1 * (i-1),
                                  custom_column_labels = labels,font.size=2.5) + coord_cartesian(clip="off")
  sub_tree_plot = annotated_basic_tree
}
sub_tree_plot


basic_tree_bars <- sub_tree_plot +
  theme(legend.position = "top",legend.text = element_text(size=text_size-0.5),
        legend.title = element_blank()) + guides(fill = guide_legend(nrow = 1))

basic_tree_bars
label_segs <- data.frame(xstart=seq(from=25+0.05,to=35-0.1,by=26/13),xend=seq(from=25+0.2,to=35-0.1,by=26/13)+1)

#adding small bars and tape labels
for(i in 1:13) {
basic_tree_bars <- basic_tree_bars + geom_segment(x=label_segs$xstart[i],y=-5,
                                                  xend=label_segs$xend[i],yend=-5,color="black")
}
basic_tree_bars

#combined
basic_tree_bars_axis
ggsave(paste0(output_dir,"MCC_subtree_with_alignment_annotated.png"), basic_tree_bars, width = 14.28, height = 15, units = "cm", dpi = 300)
svg(filename = paste0(output_dir,"MCC_subtree_with_alignment_annotated.svg"), width = 14.28, height = 15)
basic_tree_bars
dev.off()
