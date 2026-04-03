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
library(tidytree)
library(treeio)
library(phytools)

text_size = 10
output_dir = "plots/"
## ---------------------------

swap_integer_for_edit = function(integer, insert_to_integer_map){

  insert = insert_to_integer_map[insert_to_integer_map$integer == integer, "insert"]
  return(insert)
}

## ---------------------------

## define input files

tree_file = "inference_output/3-CCD_clockPerTarget_sampling_DataSet1_3000000.tree"
alignment = "pre_processed_data/1-alignment_seed1.txt"
cell_ids_file = "pre_processed_data/1-cell_ids_seed1.txt"
edit_file = "pre_processed_data/1-edit_table_sample_1.csv"
edit_to_integer_map = "pre_processed_data/1-insert_to_integer_map.csv"


# Hack from github due to package error: could not find function "offspring.tbl_tree_item"
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

## load input
tree = treeio::read.beast(file = tree_file)
dat = tree@data

stem_length = 25 - max(nodeHeights(tree@phylo))

tree_edges = tree@phylo$edge

for(node_number in 1:1999) { 
  
if((length(caper::clade.members(node_number,tree@phylo)) < 110) && (length(caper::clade.members(node_number,tree@phylo)) > 100) ) { 
  print(node_number)
  break
}
}

# 1. Generate the base plot
p <- ggtree(tree, root.position = stem_length, size=0.05) +
  geom_hilight(node=1309, fill="steelblue", alpha=0.5) +
  geom_rootedge(rootedge = stem_length, size=0.05) + 
  theme_tree2()

# 2. Extract plot data to find the visual range
d <- p$data
# Get all offspring of node 1309
clade_nodes <- tidytree::offspring(tree, 1309)

# 3. Find the tips with the extreme Y-coordinates within that clade
clade_tips <- d[d$node %in% clade_nodes & d$isTip, ]
taxa1_node <- clade_tips$node[which.min(clade_tips$y)]
taxa2_node <- clade_tips$node[which.max(clade_tips$y)]

# 4. Add the bracket using these geometrically-correct nodes
p + geom_strip(
  taxa1 = taxa1_node, 
  taxa2 = taxa2_node, 
  color = "steelblue", 
  barsize = 1,
  offset = 0.5,     # Adjust to move the bracket away from the tips
  extend = 0.2      # Adjust to make the bracket slightly taller/shorter
) +
  theme(
    axis.text.x  = element_text(size = 6),
    axis.title.x  = element_text(size = text_size),
    plot.margin = margin(5, 40, 5, 5) # Added right margin for the bracket
  ) +
  xlab("Time [d]")


# 1. Generate the base plot first
p_circular <- ggtree(tree, layout = "circular", root.position = stem_length, size = 0.05) +
  geom_hilight(node = 1309, fill = "steelblue", alpha = 0.5) +
  geom_rootedge(rootedge = stem_length, size = 0.05) + 
  theme_tree() + theme(
    plot.margin = margin(0, 0, 0, 0, "pt"),
    panel.spacing = unit(0, "pt")
  )
  

# 1. Get the list of all descendant node IDs (this returns a vector)
offspring_nodes <- tidytree::offspring(tree, 1309)

# 2. Filter the plot data 'd' to find only the TIPS within that list
# We use %in% to match the IDs
clade_tips <- d[d$node %in% offspring_nodes & d$isTip, ]

# 3. Identify the geometric start and end based on the plot angles
# It's safer to use 'node' indices for geom_strip to avoid label matching errors
taxa1_node <- clade_tips$node[which.min(clade_tips$angle)]
taxa2_node <- clade_tips$node[which.max(clade_tips$angle)]

# 4. Add the strip
p_circular <- p_circular + geom_strip(
  taxa1 = taxa1_node, 
  taxa2 = taxa2_node, 
  color = "steelblue",
  barsize = 0.8,
  offset = 0.1,
  extend = 0.2
) 


# zoom in on clade and show alignment and time uncertainty
sub =  tree_subset(tree = tree, node = 1309,  group_node = T, root_edge = T, levels_back = F)
subdat = sub@data
subtree = get.tree(sub)

subtree_stem_length = 25 - max(nodeHeights(subtree))

# sub tree with negative branch lengths
sub_tree_plot = ggtree(sub, root.position = subtree_stem_length) + theme_tree2() + geom_rootedge(rootedge = subtree_stem_length) +  
  geom_range('height_0.95_HPD', color='grey', size=3, alpha=.4)+
  geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.1, size=3) 
sub_tree_plot

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

sub_tree_plot = ggtree(sub, root.position = subtree_stem_length) + theme_tree2() + geom_rootedge(rootedge = subtree_stem_length) +  
  geom_range('height_0.95_HPD', color='grey', size=3, alpha=.4)+
  geom_nodepoint(aes(size=posterior))+ scale_size(range = c(0.1, 2), name = "Posterior support") 
  

sub_tree_plot

for (i in 1:length(targetBCs)){
  print(i)
  labels <- c(targetBCs[i],rep("",ncol(targetbc_edits_list[[i]])-1))
  annotated_basic_tree = gheatmap(sub_tree_plot, targetbc_edits_list[[i]],
                                  hjust=0,colnames_offset_y=-2,colnames_offset_x=0.2,
                                  colnames_position = "bottom",width = 0.04, offset = 1.1 * (i-1),
                                  custom_column_labels = labels,font.size=2.5,colnames_angle=-45) + coord_cartesian(clip="off") 
  sub_tree_plot = annotated_basic_tree
}
sub_tree_plot 


basic_tree_bars <- sub_tree_plot +
  theme(legend.position = "top",legend.text = element_text(size=text_size-2),
         legend.spacing.x = unit(0.01, 'cm'),legend.key.size = unit(0.2, 'cm')) + guides(fill = guide_legend(nrow = 2)) 
#basic_tree_bars
label_segs <- data.frame(xstart=seq(from=25.1,to=39.5,by= 14.5/13),xend=seq(from=25.1,to=39.5,by= 14.5/13) +0.7)

#adding small bars and tape labels
for(i in 1:13) {
basic_tree_bars <- basic_tree_bars + geom_segment(x=label_segs$xstart[i],y=-1,
                                                  xend=label_segs$xend[i],yend=-1,color="black")
}
final_tree <- basic_tree_bars  + ggtree::vexpand(.1, -1)


ggsave(paste0(output_dir,"CCD_subtree_with_alignment_annotated_sampling.pdf"), final_tree, width = 14.28, height = 15, units = "cm", dpi = 300)
