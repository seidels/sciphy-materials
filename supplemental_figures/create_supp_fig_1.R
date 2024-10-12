## ---------------------------
##
## Script name: plot_edit_outcomes_per_site
##
## Purpose of script: Plot the number of edits introduced at different sites
## based on data in https://doi.org/10.1038/s41586-022-04922-8
##
## Author: Sophie Seidel
##
## Date Created: 2022-08-22
##
## Copyright (c) Sophie Seidel, 2022
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes: data obtained from https://github.com/shendurelab/DNATickerTape/
##
##
## ---------------------------

## set working directory for Mac

library(ggplot2)
library(reshape2)

## set figure settings
figure_path = "plots/supp_fig_1.pdf"


# load the preprocessed data ---------

edit_table_by_5 = readRDS("../figure_3/data/edit_table_filtered.RDS")

edit_table_by_5$nEditedSites = rowSums(!is.na(edit_table_by_5[, c("Site1", "Site2", "Site3", "Site4", "Site5")]))

edited_sites_per_tape = ggplot(data = edit_table_by_5, aes(x=nEditedSites)) + 
  geom_bar()+
  facet_wrap(~ TargetBC)+
  labs(x = "# insertions", y = "count") +
  theme_minimal()

edited_sites_per_tape

ggsave(plot =edited_sites_per_tape, filename = figure_path)
