# Setting the initial parameter and directory
options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(ape)
library(phangorn)
library(tidyverse)

# 
TargetBC_raw = read.csv('mGASv2_Lane1_TargetBC_10X_bamExtractV2_t10.csv', 
                         stringsAsFactors = F, header = T, na.strings=c("","NA"))
TargetBC_list <- filter(TargetBC_raw, mBC != "GGGGGGGGGGGG") %>%
  filter(mBC != "AAAAAAAAAAAA") %>%
  group_by(mBC) %>%
  summarise(read_total = sum(n_reads)) %>%
  filter(read_total > 100000)
hist(TargetBC_list$read_total, xlim = c(0,1000000),breaks = 200)
nTargetBC =  nrow(TargetBC_list)


CellBC_list = unique(unlist(TargetBC_raw$cBC))
hist(TargetBC_raw$n_reads, xlim = c(0,2000),breaks = 200)

CellxTargetBC <- filter(TargetBC_raw, mBC %in% TargetBC_list$mBC) %>%
  filter(n_reads > 1000)


TargetBC_overlaps <- matrix(0, nTargetBC, nTargetBC)
colnames(TargetBC_overlaps) <- TargetBC_list$mBC
rownames(TargetBC_overlaps) <- TargetBC_list$mBC


for (CellBC in CellBC_list){
  ByCell_df <- filter(CellxTargetBC, cBC == CellBC)
  TBC_list <- unique(ByCell_df$mBC)
  if (length(TBC_list) > 1){
    for (ii in 1:(length(TBC_list)-1)){
      for (jj in (ii):length(TBC_list)){
        TargetBC_overlaps[TBC_list[ii],TBC_list[jj]] <- TargetBC_overlaps[TBC_list[ii],TBC_list[jj]] + 1
        TargetBC_overlaps[TBC_list[jj],TBC_list[ii]] <- TargetBC_overlaps[TBC_list[ii],TBC_list[jj]]
      }
    }
  }
}


TargetBC_DM <-1 / (TargetBC_overlaps+1)
diag(TargetBC_DM) <- 0
TargetBC_tree <- hclust(as.dist(TargetBC_DM), "complete")

color_set3 <- brewer.pal(n = 8, name = "Set1")
member <- as.data.frame(cutree(TargetBC_tree, k = 8))
member$TargetBC <- rownames(member)
colnames(member) <- c('group','TargetBC')

TargetBC_tree <- as.phylo(hclust(as.dist(TargetBC_DM),"complete"))
plot(TargetBC_tree,tip.col=color_set3[member$group])


pdf(file = "TargetBC_clustering.pdf")
plot(TargetBC_tree,tip.col=color_set3[member$group])
dev.off()

write.csv(member,file = "TargetBC_cell_assignment_v2.csv")
