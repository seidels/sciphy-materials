## ---------------------------
##
##
## Purpose of script: Visualise edit distribution in mGASv2 dataset
##
## Author: Sophie Seidel
##
## Date Created: 2024-05-20
##
## Copyright (c) Sophie Seidel, 2024
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


## load up the packages we will need:  

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
## ---------------------------

## load up our functions into memory
source("processing_scripts/useful_scripts_across_categories.R")

## ---------------------------
#figure settings
text_size = 10

## ---------------------------
# Data preprocessing

# input file
filtered_dat_file = "data/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS"
filtered_dat = readRDS(filtered_dat_file)

# sort targetBCs alphabetically and keep record of the ordering
targetbcs_sort = data.frame(tape_names = paste0(1:8), tape_identifier = sort(unique(filtered_dat$TargetBC)))
write.csv(x = targetbcs_sort, file = "tape_seq_identifier_map.csv", quote = F, row.names = F)

# reshape dataframe for plotting
edits_melted = melt(filtered_dat, id.vars = c("TargetBC", "Cell"), variable.name = "Sites")
edits_melted = edits_melted[edits_melted$value != "None", ]
head(edits_melted)

# Add new target bc/ tape names
edits_melted <- edits_melted %>%
  mutate(TargetBC = factor(TargetBC, levels = sort(unique(TargetBC)))) %>%
  mutate(TargetBC_new = paste0(as.numeric(TargetBC)))


edits_melted <- edits_melted %>%
  mutate(SiteNum = as.numeric(gsub("Site", "", Sites)))

## ---------------------------
##### Ask interesting questions ;P

## What is the number of edited sites? (Larger number of edited sites per targetbc --> more info)
head(filtered_dat)

g =
  ggplot(edits_melted, 
         aes(x=SiteNum, col=TargetBC_new, group=TargetBC_new)) +
  geom_point(stat = "count")+
  geom_line(stat = "count")+
  #geom_bar()+
  scale_color_brewer(palette = "Dark2") +
  theme_minimal()+ 
  ylab("Insert count") +
  xlab("Sites")+
  #coord_cartesian(ylim = c(0, 100))+
  theme(text = element_text(size = text_size), 
        legend.position = "right", legend.title = element_blank()) #axis.text.x = element_blank()

g

ggsave( filename = "plots/insert_count_per_site.png", plot = g, width = 6.8, height = 7.2, units = "cm", dpi = 300)

## ---------------------------
## What is the number of unique sequences? (Larger number of unique sequences --> more info)


# Reduce to cells with different sequences
filtered_dat_unique_cells = unique(filtered_dat[, 2:7])

nr_unique_sequences_per_targetbc = sapply(targetbcs_sort$tape_identifier, function(x){
  nrow(filtered_dat_unique_cells[which(filtered_dat_unique_cells$TargetBC ==x), ])
})

unique_seqs_per_targetbc = data.frame(targetbc = targetbcs_sort$tape_identifier, 
                                      targetbc_new = targetbcs_sort$tape_names,
                                      unique_seqs = nr_unique_sequences_per_targetbc)

g = ggplot(unique_seqs_per_targetbc, aes(x=targetbc_new, y=unique_seqs, fill=targetbc_new)) + 
  geom_col()+
  ylab("Unique seq count")+
  xlab("Tapes")+
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal()+
  theme(legend.position = "none", 
        axis.title.y = element_text(hjust = 0.2))
g
ggsave( filename = "nr_unique_seqs_per_tape.png", plot = g, width = 7.5, height = 5, units = "cm", dpi = 300)


### reshape data again
edits_melted = melt(filtered_dat_unique_cells, id.vars = c("TargetBC"))
edits_melted = edits_melted[edits_melted$value != "None", ]

# Add new target bc/ tape names
edits_melted <- edits_melted %>%
  mutate(TargetBC = factor(TargetBC, levels = sort(unique(TargetBC)))) %>%
  mutate(TargetBC_new = paste0("tape-", as.numeric(TargetBC)))

## only keep 3 targetbcs 
keep_targetbcs = c("AACTGAATGTTT", "TAACGAAGATTT" , "TTGATAACGTGA" )
edits_melted = edits_melted[which(edits_melted$TargetBC %in% keep_targetbcs), ]

## rename targetBCs
new_names = paste0("Target-BC-"
              , c(1,6, 8))
edits_melted$new_name_digit = sapply(edits_melted$TargetBC, function(x){ 
  grep(x, keep_targetbcs, value = FALSE, fixed = T)
})
edits_melted$new_name = new_names[edits_melted$new_name_digit]

# sort target bcs alphabetically
edits_melted$TargetBC = factor(edits_melted$TargetBC, levels = sort(unique(filtered_dat_unique_cells$TargetBC)))

g =
  ggplot(data = subset(x = edits_melted, 
                       edits_melted$variable %in% c("Site1", "Site2", "Site3", "Site4", "Site5")), 
         aes(x=variable)) +
  facet_grid( ~new_name )+
  
  geom_bar()+
  theme_minimal()+ 
  ylab("Insert count") +
  xlab("Insert")+
  #coord_cartesian(ylim = c(0, 100))+
  theme(axis.text.x = element_blank(), axis.title = element_text(size = text_size*2))

g =
  ggplot(edits_melted, 
         aes(x=variable)) +
  facet_grid(new_name~.) +
  geom_bar()+
  theme_minimal()+ 
  ylab("Insert count") +
  xlab("Site")+
  #coord_cartesian(ylim = c(0, 100))+
  theme(text = element_text(size = text_size*2)) #axis.text.x = element_blank()
g

ggsave( filename = "insert_count.png", plot = g, width = 4.76, height = 5, units = "cm", dpi = 300)

svg(filename = "insert_count.svg", width = 4.76, height = 5, pointsize = 24)
g
dev.off()
  