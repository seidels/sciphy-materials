
library(tracerer)
library(HDInterval)
library(ggplot2)
library(stringr)
library(scales)

## --------------------------
setwd("/Users/azwaans/Documents/Projects/TYPEWRITER/runtime_per_tape_length")

#ESS obtained for 1hr runs on the cluster for a tree of 21 tips (seed=1). ESS value found with Tracer
all_stats <- data.frame(num_target=c(5,10,20,30,40,50,60,70,80,90), ess_value_likelihood=c(16791,10160,3864,3625,1490,538,584,147,211,141))
#time in hrs
time <- 1
#rescale the runtime to obtain ESS=200
all_stats$runtime_till_ESS_200 <- time*200/all_stats$ess_value_likelihood

runtime_figure <- ggplot(all_stats, aes(x=num_target, y=runtime_till_ESS_200)) +
  geom_point() + geom_smooth(method="lm") +
  xlab("Sites per tape")+
  ylab("Total runtime until convergence [hrs]") + theme_classic()+ 
  theme(legend.position = "None",text=element_text(size = 10),axis.text.x = element_text(size=7,angle = 90),axis.text.y = element_text(size=7)) +
  scale_x_continuous(breaks = seq(0,100,by=10)) 

runtime_figure
ggsave("runtime_till_ESS_200_per_tape_length_small.png",runtime_figure, width = 7.14, height = 7.14, units = "cm", dpi = 300)
ggsave("runtime_till_ESS_200_per_tape_length_small.pdf",runtime_figure, width = 7.14, height = 7.14, units = "cm", dpi = 300)


#ESS obtained for runs on the cluster for a tree of 123 tips (seed=7). ITe total runtime is found with command stat. ESS value found with Tracer
all_stats <- data.frame(num_target=c(5,10,20,30,40,50,60,70,80),ess_value_likelihood=c(1609,1208,1411,994,2054,1327,1110,672,1109))
#time in hrs
time <-c(1.45,4.33,14.5,40,40,40,40,40,40)
#rescale the runtime to obtain ESS=200
all_stats$runtime_till_ESS_200 <- time*200/all_stats$ess_value_likelihood

runtime_figure <- ggplot(all_stats, aes(x=num_target, y=runtime_till_ESS_200)) +
  geom_point() + geom_smooth(method="lm") +
  xlab("Sites per tape")+
  ylab("Total runtime until convergence [hrs]") + theme_classic()+ 
  theme(legend.position = "None",text=element_text(size = 10),axis.text.x = element_text(size=7,angle = 90),axis.text.y = element_text(size=7)) +
  scale_x_continuous(breaks = seq(0,100,by=10)) 
  
runtime_figure

write.csv(all_stats,"supp_fig_2.csv")
ggsave("runtime_till_ESS_200_per_tape_length_medium.png",runtime_figure, width = 7.14, height = 7.14, units = "cm", dpi = 300)
ggsave("runtime_till_ESS_200_per_tape_length_medium.pdf",runtime_figure, width = 7.14, height = 7.14, units = "cm", dpi = 300)


#ESS obtained for runs on the cluster for a tree of 556 tips (seed=2). ITe total runtime is found with command stat. ESS value found with Tracer
all_stats <- data.frame(num_target=c(5,10,20,30,40,50,60,70,80),ess_value_likelihood=c(137,84,98,56,104,48,23,14,12))
#time in hrs
time <-c(16.5,52.16,75.16,75.16,75.16,75.16,75.16,75.16,75.16)
#rescale the runtime to obtain ESS=200
all_stats$runtime_till_ESS_200 <- time*200/all_stats$ess_value_likelihood

runtime_figure <- ggplot(all_stats, aes(x=num_target, y=runtime_till_ESS_200)) +
  geom_point() + geom_smooth(method="lm") +
  xlab("Sites per tape")+
  ylab("Total runtime until convergence [hrs]") + theme_classic()+ 
  theme(legend.position = "None",text=element_text(size = 10),axis.text.x = element_text(size=7,angle = 90),axis.text.y = element_text(size=7)) +
  scale_x_continuous(breaks = seq(0,100,by=10)) 

runtime_figure
ggsave("runtime_till_ESS_200_per_tape_length_large.png",runtime_figure, width = 7.14, height = 7.14, units = "cm", dpi = 300)
ggsave("runtime_till_ESS_200_per_tape_length_large.pdf",runtime_figure, width = 7.14, height = 7.14, units = "cm", dpi = 300)





#