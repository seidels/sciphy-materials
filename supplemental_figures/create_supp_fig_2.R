## ---------------------------
##
##
## Purpose of script: Plot inference results for the analysis where only 2 sites were included
##
## Author: Sophie Seidel
##
## Date Created: 2024-03-18
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


## set working directory where the log files are

data_dir = "supp_fig_2_data/"

## figure settings
text_size = 11
figure_path = "./plots/supp_fig_2.pdf"

## load up the packages we will need:

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)

#load the combined log file

log_file_2_sites <- "analysis_2sites_DataSet1.combined_b30.log"
dat_2_sites <- read.table(log_file_2_sites, header = T)

log_file_all_sites = "../figure_3/inference_output/1-combined.log"
dat_all_sites = read.table(log_file_all_sites, header = T)


# --------------------
#  plot the clock rates
# --------------------

#extract the clock rate from the tract
clock_rate_all_sites = dat_all_sites[, startsWith(x = colnames(dat_all_sites), prefix = "clock")]
clock_rate_2_sites <- dat_2_sites[, startsWith(x = colnames(dat_2_sites), prefix = "clock")]

#rename by the targetBC
names(clock_rate_2_sites) = names(clock_rate_all_sites) <- c("ATGGTAAG","ATTTATAT",
                       "ATTTGGTT", "GCAGGGTG",
                       "GTAAAGAT", "TAGATTTT",
                       "TGCGATTT", "TGGACGAC",
                       "TGGTTTTG", "TTAGATTG",
                       "TTGAGGTG","TTTCGTGA",
                       "TTCACGTA")


#reorder by median
ord_all = order(-sapply(clock_rate_all_sites, median))
clock_rate_all_sites <- clock_rate_all_sites[ord_all]

ord_2 = order(-sapply(clock_rate_2_sites, median))
clock_rate_2_sites <- clock_rate_2_sites[ord_2]

ordered_names_all <- names(clock_rate_all_sites)
ordered_names_2 <- names(clock_rate_2_sites)

#add a prior column
#clock_rate <- cbind(clock_rate,Prior=rlnorm(nrow(clock_rate), meanlog = -2, sdlog = 0.5))
clock_rate_all_sites_long <- pivot_longer(clock_rate_all_sites, seq(1,ncol(clock_rate_all_sites)))
clock_rate_all_sites_long$data = "all sites"
clock_rate_2_sites_long <- pivot_longer(clock_rate_2_sites, seq(1,ncol(clock_rate_2_sites)))
clock_rate_2_sites_long$data = "2 sites"

#order columns
clock_rate_all_sites_long <- mutate(clock_rate_all_sites_long,name = fct_relevel(name,ordered_names_all))
clock_rate_2_sites_long <- mutate(clock_rate_2_sites_long, name = fct_relevel(name, ordered_names_2))

clock = rbind(clock_rate_all_sites_long, clock_rate_2_sites_long)


p_clock_pos <- ggplot(clock_rate_all_sites_long, aes(x=name,value,fill=name)) +
  facet_grid(data~.)+
  theme_bw() +
  geom_violin(draw_quantiles =  c(0.5)) +
  xlab("Tape") +
  ylab(expression("Posterior edit rate [" * d^-1 * "]"))+

  theme(legend.position = "none", axis.title = element_text(size = text_size)) +
  scale_fill_manual(values=c(rep("#5CA17D",13),"#E1E1F7")) +
  coord_cartesian(
    ,expand=FALSE) + theme(axis.title = element_text(size= text_size *2),
                           axis.text = element_text(size = (text_size -2)*2),
                           panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                           strip.background  = element_rect(colour="black", fill="white"),
                           panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 90) )

p_clock_pos
ggsave("compare_clock_all_sites.png", p_clock_pos, width = 14.28, height = 9, units = "cm", dpi = 300)
svg("compare_clock_all_sites.svg", width = 14.28, height = 9)
p_clock_pos
dev.off()

p_clock_pos <- ggplot(clock_rate_2_sites_long, aes(x=name,value,fill=name)) +
  facet_grid(data~.)+
  theme_bw() +
  geom_violin(draw_quantiles =  c(0.5)) +
  xlab("Tape") +
  ylab(expression("Posterior edit rate [" * d^-1 * "]"))+

  theme(legend.position = "none", axis.title = element_text(size = text_size)) +
  scale_fill_manual(values=c(rep("#5CA17D",13),"#E1E1F7")) +
  coord_cartesian(
    ,expand=FALSE) + theme(axis.title = element_text(size= (text_size)*2),
                           axis.text = element_text(size = (text_size -2)*2),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           strip.background  = element_rect(colour="black", fill="white"),
                           panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 90) )

p_clock_pos

svg("compare_clock_2_sites.svg", width = 14.28, height = 9, pointsize = )
p_clock_pos
dev.off()
ggsave("compare_clock_2_sites.png", p_clock_pos, width = 14.28, height = 9, units = "cm", dpi = 300)



# Test the average of the medians

## Among all sites
tapes_4x = ordered_names_all[2:4]
tapes_5x = ordered_names_all[5:length(ordered_names_all)]

medians_tapes4x_allsites = apply(X = clock_rate_all_sites[, tapes_4x], MARGIN = 2, FUN = median)
mean_tapes4x_allsites = mean(medians_tapes4x_allsites)
mean_tapes4x_allsites
sd(medians_tapes4x_allsites)

medians_tapes5x_allsites = apply(X = clock_rate_all_sites[, tapes_5x], MARGIN = 2, FUN = median)
mean_tapes5x_allsites = mean(medians_tapes5x_allsites)
mean_tapes5x_allsites
sd(medians_tapes5x_allsites)

## Among 2 sites
medians_tapes4x_2sites = apply(X = clock_rate_2_sites[, tapes_4x], MARGIN = 2, FUN = median)
mean_tapes4x_2sites = mean(medians_tapes4x_2sites)
mean_tapes4x_2sites
sd(medians_tapes4x_2sites)

medians_tapes5x_2sites = apply(X = clock_rate_2_sites[, tapes_5x], MARGIN = 2, FUN = median)
mean_tapes5x_2sites = mean(medians_tapes5x_2sites)
mean_tapes5x_2sites
sd(medians_tapes5x_2sites)




