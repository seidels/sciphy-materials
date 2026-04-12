# libs
library(ggplot2)
library(tidyverse)

figure_dir = "site-simulations/plots/"

make_clockrate_long <- function(dat) {
  clock_levels <- dat %>%
    select(starts_with("clockRate_")) %>%
    names() 
  
  dat %>%
    select(all_of(clock_levels)) %>%
    pivot_longer(cols = everything()) %>%
    mutate(
      name = factor(name, levels = clock_levels),
      tape_4x = name %in% c("clockRate_10", "clockRate_11", "clockRate_12")
    )
}

# Vector of seed values
seeds <- c(1, 4,5)

# Construct paths and process 
dat_norm <- map_dfr(seeds, function(seed) {
  path <- paste0("supp_fig_7_ratediff/ratediff0.5/alignments_seed_", seed, "/infer_normal.100000.combined.log")
  dat <- read.delim(path)
  dat_long <- make_clockrate_long(dat)
  dat_long$rate_diff <- "No rate diff."
  
  if(seed == 1 ) { 
    dat_long$seed <- "Tree 1"}
  if(seed == 4 ) { 
    dat_long$seed <- "Tree 2"}
  if(seed == 5 ) { 
    dat_long$seed <- "Tree 3"}
  dat_long
})

dat_diff <- map_dfr(seeds, function(seed) {
  path <- paste0("supp_fig_7_ratediff/ratediff0.2/alignments_seed_", seed, "/infer_ratediff.100000.combined.log")
  dat <- read.delim(path)
  dat_long <- make_clockrate_long(dat)
  dat_long$rate_diff <- "Rate diff."
  
  if(seed == 1 ) { 
    dat_long$seed <- "Tree 1"}
  if(seed == 4 ) { 
    dat_long$seed <- "Tree 2"}
  if(seed == 5 ) { 
    dat_long$seed <- "Tree 3"}
  dat_long
})
dat_all = rbind(dat_norm, dat_diff)

# colour top 3 median clock rates
top3_names <- dat_all %>%
  group_by(seed, rate_diff, name) %>%
  summarise(median_value = median(value), .groups = "drop") %>%
  arrange(seed, rate_diff, desc(median_value)) %>%
  group_by(seed, rate_diff) %>%
  slice_head(n = 3) %>%
  ungroup()

#Add a new variable to dat_all marking if a row is among top 3
dat_all <- dat_all %>%
  left_join(top3_names %>% select(seed, rate_diff, name) %>% mutate(top3 = TRUE),
            by = c("seed", "rate_diff", "name")) %>%
  mutate(top3 = ifelse(is.na(top3), FALSE, top3))


p <-ggplot(dat_all, aes(x=name, y=value, colour = tape_4x, fill=top3)) + 
  geom_violin()+
  facet_grid(rate_diff~ seed) +
  scale_x_discrete(labels = c("clockRate_1" = "1", 
                              "clockRate_2" = "2",
                              "clockRate_3" = "3",
                              "clockRate_4" = "4",
                              "clockRate_5" = "5",
                              "clockRate_6" = "6",
                              "clockRate_7" = "7",
                              "clockRate_8" = "8",
                              "clockRate_9" = "9",
                              "clockRate_10" = "10",
                              "clockRate_11" = "11",
                              "clockRate_12" = "12"
                              )
  )+
  geom_hline(yintercept = 0.1, colour = "darkgreen")+
  scale_color_manual(values = c("black", "#0072B2"))+
  scale_fill_manual(values = c("white", "#E69F00"))+
  theme_bw() +
  theme(legend.position = "top", panel.grid.major.y = element_line(color = "#ebebeb", linewidth = 0.2), 
        panel.grid.minor.y = element_line(color = "#ebebeb", linewidth = 0.2), 
        panel.grid.major.x = element_blank(),
    panel.background = element_blank(),
     strip.text = element_text(size = 10),
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 0.5)
  ) +  
  
  labs(y = expression("Posterior editing rate [" * d^-1 * "]"),x="Tape index") 

ggsave("plots/supp_fig_7.pdf", plot = p,height = 14.28, units = "cm", dpi = 300)
