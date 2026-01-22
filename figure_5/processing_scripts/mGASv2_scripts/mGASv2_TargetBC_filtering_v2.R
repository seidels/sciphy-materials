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

# Reading in the cell-by-TAPE
Cell_By_TAPE = read.csv('mGASv2_Lane1_CellByTape_10X_bamExtractV2_t3_collapse.csv', 
                        stringsAsFactors = F, header = T, na.strings=c("","NA"))
Cell_By_TAPE$Cell <- substr(Cell_By_TAPE$Cell,1,16)

cell_list_TAPE_total <- unique(Cell_By_TAPE$Cell)
mGASv2_group1_TargetBC <- filter(member, group == 2)
mGASv2_Cell_By_TAPE = filter(Cell_By_TAPE, TargetBC %in% mGASv2_group1_TargetBC$TargetBC)

# In EBv2, lane1, group1, there are 31 TargetBC. Counting # of cells with TargetBC 
mGASv2_TBC_count <- mGASv2_Cell_By_TAPE
mGASv2_TBC_count$count <- 1
mGASv2_TBC_count_grouped <- mGASv2_TBC_count %>%
  group_by(TargetBC)%>%
  summarise(
    TBC_count = sum(count)
  )



# Calculating entropy associated with each TargetBC

TBC_all_entropy <- data.frame()
TBC_samples_all <- data.frame()
for (TBC in mGASv2_TBC_count_grouped$TargetBC){
  TBC_temp <- filter(mGASv2_Cell_By_TAPE, TargetBC == TBC) %>%
    mutate(Sites = paste(Site1,Site2,Site3,Site4,Site5,Site6))
  TBC_freq <- table(TBC_temp$Sites)/dim(TBC_temp)[1]
  TBC_entropy <- (-1)*sum(TBC_freq*log(TBC_freq))
  #TBC_all_entropy[TBC] <- TBC_entropy
  TBC_samples <- arrange(data.frame(table(TBC_temp$Sites)), desc(Freq))[1:10,]
  TBC_samples$TargetBC <- TBC
  TBC_samples$Entropy <- TBC_entropy
  TBC_all_entropy <- bind_rows(TBC_all_entropy, data.frame(TBC, TBC_entropy))
  TBC_samples_all <- bind_rows(TBC_samples_all, TBC_samples)
}
colnames(TBC_all_entropy) <- c('TargetBC','Entropy')

# Filtering by entropy and recovery counts
mGASv2_TBC_count_grouped <- left_join(mGASv2_TBC_count_grouped, TBC_all_entropy, by = 'TargetBC')
