
library(tracerer)
library(HDInterval)
library(ggplot2)
library(stringr)
library(scales)

## --------------------------

all_stats <- read.csv("supp_fig_runtime_data/supp_fig_2.csv")

runtime_figure <- ggplot(all_stats, aes(x=num_target, y=runtime_till_ESS_200)) +
  geom_point() + geom_smooth(method="lm") +
  xlab("Sites per tape")+
  ylab("Total runtime until convergence [hrs]") + theme_classic()+ 
  theme(legend.position = "None",text=element_text(size = 10),axis.text.x = element_text(size=7,angle = 90),axis.text.y = element_text(size=7)) +
  scale_x_continuous(breaks = seq(0,100,by=10)) 

runtime_figure

ggsave("plots/supp_fig_2.pdf",runtime_figure, width = 14.28, height = 14.28, units = "cm", dpi = 300)
