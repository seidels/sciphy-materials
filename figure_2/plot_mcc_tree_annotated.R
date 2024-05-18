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

## set working directory for Mac

setwd("~/Projects/typewriter_analysis")    # personal working directory

## load up the packages we will need:  (uncomment as required)
library(ggtree)
library(ggplot2)
library(treeio)

text_size = 12
## ---------------------------

swap_integer_for_edit = function(integer, insert_to_integer_map){

  insert = insert_to_integer_map[insert_to_integer_map$integer == integer, "insert"]
  return(insert)
}

## ---------------------------

# output dir
output_dir = "paper_figures/figure_2/"

## define input files

tree_file = "results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/MCC_on_thinned4000000trees.tree"
alignment = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/alignment_seed1.txt"
cell_ids_file = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/cell_ids_seed1.txt"
edit_file = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/edit_table_sample_1.csv"
edit_to_integer_map = "data/cell_culture/insert_to_integer_map.csv"

## load input
tree = read.beast(file = tree_file)
insert_to_integer_map = read.csv(edit_to_integer_map)

cell_ids = read.csv(cell_ids_file, header = F)
cell_ids$tip_label = 0:999 # keep numbering of tips as used in alignment
cell_ids = cell_ids[order(cell_ids$V1), ]


edits = read.csv(edit_file)

## create targetbc matrices
ctr = 1
targetBCs = unique(edits$TargetBC)
targetbc_edits_list =  vector(mode = "list", length = length(targetBCs))

for (targetBC in targetBCs){

  print(targetBC)
  edits_tbc = edits[ edits$TargetBC == targetBC, ]

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
basic_tree = ggtree(tree,layout="ellipse")


for (i in 1:length(targetBCs)){
  print(i)
  labels <- c(targetBCs[i],rep("",ncol(targetbc_edits_list[[i]])-1))
  annotated_basic_tree = gheatmap(basic_tree, targetbc_edits_list[[i]],
                                  hjust=0,colnames_offset_y=-10,colnames_offset_x=0,
                                  colnames_position = "bottom",width = 0.07, offset = 2.0 * (i-1),
                                  custom_column_labels = labels,font.size=2.5) + coord_cartesian(clip="off")
  basic_tree = annotated_basic_tree
}

basic_tree_bars <- basic_tree +
  theme(legend.position = "top",legend.text = element_text(size=text_size),
        legend.title = element_blank()) + guides(fill = guide_legend(nrow = 1))

label_segs <- data.frame(xstart=seq(from=25+0.05,to=51-0.1,by=26/13),xend=seq(from=25+0.2,to=51-0.1,by=26/13)+1)

#adding small bars and tape labels
for(i in 1:13) {
basic_tree_bars <- basic_tree_bars + geom_segment(x=label_segs$xstart[i],y=-5,xend=label_segs$xend[i],yend=-5,color="black")
}
#adding a makeshift axis
basic_tree_bars_axis <- basic_tree_bars +
  geom_segment(x=0,y=-5,xend=24.8,yend=-5,color="black") +
  geom_text(x=0,y=-15,label="0",size=8) +
  geom_text(x=5,y=-15,label="5",size=8)  +
  geom_text(x=10,y=-15,label="10",size=8)  +
  geom_text(x=15,y=-15,label="15",size=8)  +
  geom_text(x=20,y=-15,label="20",size=8) +
  coord_cartesian(ylim=c(-20,1000),clip="off")

#combined
ggsave(paste0(output_dir,"MCC_with_alignment_annotated.png"), basic_tree_bars_axis, width = 14.28, height = 9.5, units = "cm", dpi = 300)

