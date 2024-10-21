## ---------------------------
##
## Script name: create_supp_figs_7_8
##
## Purpose of script: Implements a simple R simulator of SciPhy data. Here, used to visualize the distributions of edits that can be expected assuming Poisson distributed number of edits along different population scenarios
##
## Author: Antoine Zwaans
##
## Date Created: 2022-09-13
##
## Copyright (c) Antoine Zwaans, 2022
## Email: antoine.zwaans@bsse.ethz.ch
##


## load up the packages we will need:  
library(TreeSim)
library(ggtree)
library(cowplot)

######################################
## helper functions for the simulation
######################################

#lambda of 0.2 allows to have close to ~ 5 edits along the tree
#draw edits numbers along each tree edge == cells, along a tree, assuming a rate lambda
#draw inserts types according to the insert frequencies
draw_edits_along_edges <- function(tree,lambda=0.20,insert_types,insert_freq) {
  num_edit_edges <- sapply(tree$edge.length, function(x) {rpois(1,x*lambda)})
  edit_types_edges <- sapply(num_edit_edges, function(x) {sample(insert_types, x, replace = TRUE, prob = insert_freq)})
  edit_types_edges <- sapply(edit_types_edges, function(x) if(identical(x, character(0))) "" else x)
  return_list <- list("num_edit"=num_edit_edges,"edits"=edit_types_edges)
  return(return_list)
}

#Get the total number of edits obtained per cell sampled at "present". 
#For each tip we sum the number of edits drawn onedge edge leading to the tip
get_edits_per_sampled_cell <- function(tree,num_edit_edges) {
  total_branch_edits <- rep(0,tree$Nnode)
  for(i in 1:length(tree$tip)) {
    start_edge <- i
    sum_path <- 0
    while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
      sum_path <- sum_path + num_edit_edges[tree$edge[,2] == start_edge]
      start_edge = tree$edge[tree$edge[,2] == start_edge,1]
      
    }
    #Specifying a target length of 5 here. Replace for any other length.  
    total_branch_edits[i] <- min(5,sum_path)
    
  }
  return(data.frame(edits=total_branch_edits))
}

##Get full barcode alignment at tips based on the ones drawn on each edge. 
## traverse the tree based on edge numbers
get_sequence_per_sampled_cell <- function(tree,edit_types_edges) {
  full_sequences <- c()
  for(i in 1:length(tree$tip)) {
    start_edge <- i
    sequence_path <- ""
    while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
      sequence_path <- paste(c(sequence_path,unlist(edit_types_edges[tree$edge[,2] == start_edge])),collapse="")
      start_edge = tree$edge[tree$edge[,2] == start_edge,1]
      
    }
    #this is to enforce a tape length of at most 5. Here we insert trinucleotides, so restricting to at most 15 nucleotides total
    print(sequence_path)
    full_sequences[i] <- min(5,nchar(sequence_path)/3)
    
    
  }
  return(full_sequences)
}

# find tree root age 
tree_height_calc <- function(tree) { 
  start_edge <- 1
  sum_path <- 0
  while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
    sum_path <- sum_path + tree$edge.length[tree$edge[,2] == start_edge]
    start_edge = tree$edge[tree$edge[,2] == start_edge,1]
    
  }
  return(sum_path)
  
}


####################################################
#input the HEK293T dataset for ~realistic simulations
####################################################

##input format
edit_table_by_5 = read.csv("../figure_3/data/supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))

bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]))) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$Freq <- bulk_insert_count$Freq/sum(bulk_insert_count$Freq)
##this is used to draw the trinucleotides in sequence simulations



################################################
#Tree simulation: simulating trees with 500 tips. 
################################################

#tree with high branching towards the present: in [0,0.3] birth rate is 10x what it is in the rest of the population history
#for high branching in present: 
#birth: c(10,1)
#death: c(0,0.3)
#mass extinction fraction survival  c(0.6,0.1)
#times c(0,0.3)
tree_present <- sim.rateshift.taxa(n=500,numbsim=1,c(10,1),c(0,0.3),
                   c(0.6,0.1),c(0,0.3),complete=FALSE)[[1]]

#for high branching in past: add death in the first epoch with respect to scenario above.  
#birth: c(10,1)
#death: c(0.3,0.3)
#mass extinction fraction survival  c(0.6,0.1)
#times c(0,0.3)
tree_past <- sim.rateshift.taxa(n=500,numbsim=1,c(1,10),c(0.3,0.3),
                            c(0.6,0.1),c(0,0.3),complete=FALSE)[[1]]

#tree with a constant birth rate
tree_even <- sim.rateshift.taxa(n=500,numbsim=1,c(10,10),c(0,0.3),
                                c(0.6,0.1),c(0,0.3),complete=FALSE)[[1]]

#find tree heights
present_height <- tree_height_calc(tree_present)
past_height <- tree_height_calc(tree_past)


#scale the tree heights in tree_past to tree_present for comparability
tree_past$edge.length <- tree_past$edge.length * (present_height/past_height)


######## SIMULATING EDIT NUMBERS AND SEQUENCES ########

### Tree with branching towards present ###

#Simulating the distribution of total number of edits at tips on the tree with high birth towards the present
#Simulating 13 barcodes with 5 positions each (default) 
total_branch_edits_present1 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2, insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present2 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present3 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present4 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present5 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present6 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present7 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present8 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present9 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present10 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present11 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present12 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present13 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)

#total at the tips is the sum for all barcodes
summed_total_present <- total_branch_edits_present1 + 
  total_branch_edits_present2 + 
  total_branch_edits_present3 + 
  total_branch_edits_present4 + 
  total_branch_edits_present5 + 
  total_branch_edits_present6 + 
  total_branch_edits_present7 + 
  total_branch_edits_present8 +
  total_branch_edits_present9 +
  total_branch_edits_present10 +
  total_branch_edits_present11 +
  total_branch_edits_present12 +
  total_branch_edits_present13

#plot the simulated number of edits/cells as a histogram + density plot
density_plot_present <- ggplot(data=summed_total_present,aes(x=edits)) + 
                        geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + 
                        xlab("Number of edits per cell") + 
                        ylab("Density") +
                        geom_density(alpha=.2, fill="#FF6666") 

#plot tree facing the histogram
plot_grid(ggtree(tree_present),density_plot_present)

#save plot
ggsave("supp_fig_7.pdf",path="plots/", width=25,height= 18, units = "cm") 
#EXAMPLE: simulate tape alignment along the tree
simulation_present <- draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)
#extract sequences and edits
sim_seqs_present <- get_sequence_per_sampled_cell(tree_present,simulation_present$edits)


### Tree with branching towards past ###

#Simulating distribution of total number of edits at tips on the tree with high branching towards the past
#Simulating 13 barcodes with 5 positions each (default) 
total_branch_edits_past1 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past2 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past3 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past4 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past5 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past6 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past7 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past8 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past9 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past10 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past11 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past12 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past13 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,lambda= 0.2,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)

#total at the tips is the sum for all barcodes
summed_total_past <- total_branch_edits_past1 + 
  total_branch_edits_past2 + 
  total_branch_edits_past3 + 
  total_branch_edits_past4 + 
  total_branch_edits_past5 + 
  total_branch_edits_past6 + 
  total_branch_edits_past7 + 
  total_branch_edits_past8 +
  total_branch_edits_past9 +
  total_branch_edits_past10 +
  total_branch_edits_past11 +
  total_branch_edits_past12 +
  total_branch_edits_past13

#plot the simulated number of edits/cells as a histogram + density plot
density_plot_past <- ggplot(data=summed_total_past,aes(x=edits)) + 
                     geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + 
                     xlab("Number of edits per cell") + 
                     ylab("Density") +
                     geom_density(alpha=.2, fill="#FF6666") 

#plot tree facing the histogram
plot_grid(ggtree(tree_past),density_plot_past)

#save plot
ggsave("supp_fig_8.pdf",path="plots/", width=25,height= 18, units = "cm")
