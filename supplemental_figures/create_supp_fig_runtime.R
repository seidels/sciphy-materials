
library(tracerer)
library(HDInterval)
library(ggplot2)
library(stringr)
library(scales)

## --------------------------
all_stat <- read.csv(file = paste0("supp_fig_runtime_data/combined_runtime_stats.csv"))
all_stat$runtime_till_ESS_200 <- all_stat$runtime * 200/(all_stat$ess_likelihood)


runtime_figure <- ggplot(all_stat, aes(x=num_tip, y=runtime_till_ESS_200/60)) +
  geom_point() + geom_smooth(size=0.5) +
  xlab("Number of cells")+
  ylab("Total runtime until convergence [hrs]") + theme_classic()+ 
  theme(legend.position = "None",text=element_text(size = 10),axis.text.x = element_text(size=7,angle = 90),axis.text.y = element_text(size=7)) +
  scale_y_continuous(trans = sqrt_trans(), 
                     breaks = c(0,2,4,6,8,10,12,14,16,18,20)^2) +
  scale_x_continuous(breaks = seq(0,650,by=50)) 
  
ggsave("supp_fig_runtime.pdf",runtime_figure, width = 7.14, height = 7.14, units = "cm", dpi = 300)
