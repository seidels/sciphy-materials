## ---------------------------
##
## Script name: plot_inference_results_annotated
##
## Purpose of script: Plot estimates from SciPhy on HEK293 cell culture data, clock per target, with color annotations matching the tree.
##
## Author: Antoine Zwaans
##
## Date Created: 2023-10-13
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory where the log files are

setwd("~/Projects/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells")
pic_dir = "~/Projects/typewriter_analysis/paper_figures/figure_2/"

# set figure settings
text_size=10
label_size = 12
## load up the packages we will need:

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)
library(magick)

#load the combined log file

typewriter_file <- "combined.log"
typewriter <- read.table(typewriter_file, header = T)

# -----------------------------
#  plot insertion probabilities
# -----------------------------

#extract insertion probabilities from the full dataframe
insert_probs <- typewriter[,c(7:25)]
trinucleotides_names <- c("CAT","CCG","GCC","ATA","GAT","ACG","ACA","TCG","TAT","GCT","CTA","TGT","AGA","TAA","CAG","TAG","GAG","ACC","GCG")

#reorder by median
names(insert_probs) <- trinucleotides_names
insert_probs <- insert_probs[order(-sapply(insert_probs, median))]
insert_probs <- bind_cols(insert_probs,Prior=rdirichlet(nrow(insert_probs), rep(1.5,19))[,2])
insert_probs_medians <- sapply(insert_probs,median)
insert_probs_low <- as.numeric(sapply(insert_probs,function(x) {quantile(x, 0.025)}))
insert_probs_up <-  as.numeric(sapply(insert_probs,function(x) {quantile(x, 0.975)}))
datafra <- data_frame(name=c(trinucleotides_names,"Prior"),median=insert_probs_medians,low=insert_probs_low,up=insert_probs_up)


datafra$name <- factor(datafra$name, levels = datafra$name)

#assign standard colors to the trinucleotides alphabetical to match the tree plot.
colors <- c(hue_pal()(19),"#E1E1F7")
names(colors) <- c(sort(trinucleotides_names),"Prior")
p_inserts <-  ggplot(datafra) +
  geom_bar(aes(x=name,y=median,fill=name),stat="identity",colour="black") +
  xlab("Insert") +
  theme_classic() +
  ylab("Posterior insert probability") + xlab("Insert") +
  theme(legend.position = "none" ) +
  scale_fill_manual(values=colors) +
  #theme(plot.title = element_text(hjust=0.5)) +
  geom_errorbar(aes(x=name,ymin=low, ymax=up), width=.2,
                position=position_dodge(.9)) +
  #ggtitle("Insert probabilities") +

  coord_cartesian(ylim=c(0,0.175),expand = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.3, vjust = 0.2,
                                   margin = margin(l = 20, r = 20), size = (text_size-4)),
        axis.text.y = element_text(size = 5, margin = margin(r = -0.1, l=0)),
        axis.title.y = element_text(size = text_size),

                                                          panel.grid.minor = element_blank(),
                                                          panel.border = element_blank(),
                                                          panel.background = element_blank(),
        panel.grid.major.x = element_blank(),axis.title.x = element_blank())


ggsave(paste0(pic_dir, "insert_probs.png"), p_inserts, width = 4.76, height = 7.14, units = "cm", dpi = 300)
svg(filename = paste0(pic_dir, "insert_probs.svg"), width = 4.76, height = 7.14, )
p_inserts
dev.off()
  # --------------------
#  plot the clock rates
# --------------------

#extract the clock rate from the tract
clock_rate <- typewriter[,c(29:41)]

#rename by the targetBC
names(clock_rate) <- c("ATGGTAAG","ATTTATAT",
                       "ATTTGGTT", "GCAGGGTG",
                       "GTAAAGAT", "TAGATTTT",
                       "TGCGATTT", "TGGACGAC",
                       "TGGTTTTG", "TTAGATTG",
                       "TTGAGGTG","TTTCGTGA",
                       "TTCACGTA")



#reorder by median
clock_rate <- clock_rate[order(-sapply(clock_rate, median))]
ordered_names <- names(clock_rate)

#add a prior column
clock_rate <- bind_cols(clock_rate,Prior=rlnorm(nrow(clock_rate), meanlog = -2, sdlog = 0.5))
clock_rate_long <- pivot_longer(clock_rate,seq(1,ncol(clock_rate)))

#order columns
clock_rate_long <- mutate(clock_rate_long,name = fct_relevel(name,ordered_names,"Prior"))

#annotating the truncated with a different colior pattern
# the 4 truncation
# manually annotate colors: darkest for most truncated
# (targetBC == "TGGACGAC") | (targetBC == "TGGTTTTG") | (targetBC == "TTTCGTGA" color #5CA17D
# truncation
# "TTCACGTA" color "#3C614F"
segments <- data.frame(xstart=c(2,5),xstop=c(4,13),ystart=c(0.25,0.21),ystop=c(0.25,0.21))
p_clock_pos <- ggplot(clock_rate_long,aes(x=name,value,fill=name)) +
  theme_classic() +
  geom_violin(draw_quantiles =  c(0.5)) +
  ylab(parse(text = paste0('"Posterior edit rate "'))) +

  theme(legend.position = "none") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3))+
  scale_fill_manual(values=c(c("#0A0A82",rep("#404085",3),rep("#7C7CA3",9)),"#E1E1F7")) +
  coord_cartesian(
    ylim=c(0.1,0.37),expand=FALSE) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = 6, margin = margin(r = 0, l=-0.5)),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = (text_size-4), margin= margin(l = -0.5, r = -0.5, b=-3)), #vjust = 0.3, hjust = 0.5),
        axis.title.x = element_blank() ) +

  annotate("text", x = 1, y = 0.35, label = "2X") +
  annotate("text", x = 3, y = 0.27, label = "4X") +
  annotate("text", x = 9, y = 0.24, label = "5X") +
  #ggtitle("Editing rate per tape")+
  #theme(plot.title = element_text(hjust=0.5)) +
  geom_segment(data=segments,aes(x=xstart,xend=xstop,y=ystart,yend=ystop),inherit.aes = FALSE)


#ggsave(paste0(pic_dir,"clock_rate.png"), p_clock_pos, width = 4.76, height = 7.14, units = "cm", dpi = 300)
ggsave(paste0(pic_dir,"clock_rate.png"), p_clock_pos, width = 9.52, height = 4.2, units = "cm", dpi = 300)
  # --------------------
#  plot the growth rate
# --------------------

#extract and substract the birth and death rates
bd_rates <- typewriter[,"birthRate"] - typewriter[,"deathRate"]

#add a prior column
bd_rates <- bind_cols(Rate=bd_rates,Prior= rlnorm(length(bd_rates), meanlog = -0.6, sdlog = 1) - rlnorm(length(bd_rates), meanlog = -2, sdlog = 1) )
bd_rates_long <- pivot_longer(bd_rates,seq(1,ncol(bd_rates)))

#order columns
bd_rates_long <- mutate(bd_rates_long,name = fct_relevel(name,"Rate","Prior"))

p_growth <- ggplot(bd_rates_long, aes(x=name,value,fill=name)) +
  theme_classic() +
  geom_violin(draw_quantiles = 0.5) +
  theme(legend.position = "none" ) +
  ylab(parse(text = paste0('"Posterior growth rate "', '(~day^-1)'))) +
  coord_cartesian(ylim = c(-0.1,0.9),expand = TRUE) +
  ggtitle("Growth rate")+
  #theme(plot.title = element_text(hjust=0.5)) +
  scale_fill_manual(values=c("white","#E1E1F7")) + theme(text=element_text(size = text_size),panel.grid.minor = element_blank(),
                                                         panel.border = element_blank(),
                                                         panel.background = element_blank(),panel.grid.major.x = element_blank(),axis.title.x = element_blank())


ggsave("growth_rate.png", p_growth, width = 15, height = 10, units = "cm", dpi = 300)


#extract and substract the birth and death rates
bd_rates <- typewriter[,"birthRate"]

#add a prior column
bd_rates <- bind_cols(Rate=bd_rates,Prior= rlnorm(length(bd_rates), meanlog = -0.6, sdlog = 1)  )
bd_rates_long <- pivot_longer(bd_rates,seq(1,ncol(bd_rates)))

#order columns
bd_rates_long <- mutate(bd_rates_long,name = fct_relevel(name,"Rate","Prior"))
p_birth <- ggplot(bd_rates_long, aes(x=name,value,fill=name)) +
  theme_bw() +
  geom_violin(draw_quantiles = 0.5) +
  theme(legend.position = "none" ) +
  xlab("Division") +
  ylab(parse(text = paste0('"Posterior division rate "', '(~day^-1)'))) +
  coord_cartesian(ylim = c(-0.1,0.9),expand = TRUE) +
  scale_fill_manual(values=c("#5CA17D","#E1E1F7")) + theme(text=element_text(size = text_size),panel.grid.minor = element_blank(),
                                                           panel.border = element_blank(),
                                                           panel.background = element_blank(),panel.grid.major.x = element_blank())


#extract and substract the birth and death rates
bd_rates <- typewriter[,"deathRate"]

#add a prior column
bd_rates <- bind_cols(Rate=bd_rates,Prior= rlnorm(length(bd_rates), meanlog = -2, sdlog = 1) )
bd_rates_long <- pivot_longer(bd_rates,seq(1,ncol(bd_rates)))

#order columns
bd_rates_long <- mutate(bd_rates_long,name = fct_relevel(name,"Rate","Prior"))
p_death <- ggplot(bd_rates_long, aes(x=name,value,fill=name)) +
  theme_bw() +
  geom_violin(draw_quantiles = 0.5) +
  theme(legend.position = "none" ) +
  xlab("Death") +
  ylab(parse(text = paste0('"Posterior death rate "', '(~day^-1)'))) +
  coord_cartesian(ylim = c(-0.1,0.9),expand = TRUE) +
  scale_fill_manual(values=c("#5CA17D","#E1E1F7")) + theme(text=element_text(size = text_size),panel.grid.minor = element_blank(),
                                                           panel.border = element_blank(),
                                                           panel.background = element_blank(),panel.grid.major.x = element_blank())



combined <- cowplot::plot_grid(p_birth+ theme(axis.line = element_line(colour = "black")),p_death+ theme(axis.line = element_line(colour = "black")),nrow = 1)



ggsave("birth_death_rates.png", combined, width = 50, height = 15, units = "cm", dpi = 300)



intro_to_dat = cowplot::draw_image(image = paste0(pic_dir, "monoclonal_expansion_HEK293T.svg"), scale = 1)
# ---------------------------
#  plot all estimates aligned
# ---------------------------
top_row <- cowplot::plot_grid(p_clock_pos+ theme(axis.line = element_line(colour = "black")), p_growth+ theme(axis.line = element_line(colour = "black")), ncol=2,align = "h",rel_widths = c(1,0.5),label_size=22)
combined <- cowplot::plot_grid(top_row,p_inserts+ theme(axis.line = element_line(colour = "black")),nrow = 2,rel_heights = c(1,0.6),label_size=text_size)
ggsave("combined_estimates_vertical.png", combined, width = 50, height = 35, units = "cm", dpi = 300)
combined <- cowplot::plot_grid(p_inserts+
                                 theme(axis.line = element_line(colour = "black")),p_clock_pos+
                                 theme(axis.line = element_line(colour = "black")), p_growth+
                                 theme(axis.line = element_line(colour = "black")),
                               nrow = 1,align = "h",
                               label_size=text_size,rel_widths = c(1,1.5,0.7))
ggsave("~/Projects/typewriter_analysis/paper_figures/figure_2/combined_estimates_horizontal.png",
       combined, width = 14.28, height = 7.14, units = "cm", dpi = 300)

row_1 <- cowplot::plot_grid(intro_to_dat + theme(axis.line = element_line(colour = "black")) +
                              p_inserts + theme(axis.line = element_line(colour = "black")),
                            p_clock_pos,
                            nrow = 1,align = "h",
                               label_size=text_size,rel_widths = c(1,1,1))
ggsave("~/Projects/typewriter_analysis/paper_figures/figure_2/row_1.png",
       row_1, width = 14.28, height = 7.14, units = "cm", dpi = 300)










