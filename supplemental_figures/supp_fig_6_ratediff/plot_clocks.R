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

# Construct paths and process all
dat_norm <- map_dfr(seeds, function(seed) {
  path <- paste0("site-simulations/xmls/ratediff0.5/alignments_seed_", seed, "/logs/infer_normal.100000.combined.log")
  dat <- read.delim(path)
  dat_long <- make_clockrate_long(dat)
  dat_long$rate_diff <- FALSE
  dat_long$seed <- seed
  dat_long
})

dat_diff <- map_dfr(seeds, function(seed) {
  path <- paste0("site-simulations/xmls/ratediff0.2/alignments_seed_", seed, "/logs/infer_ratediff.100000.combined.log")
  dat <- read.delim(path)
  dat_long <- make_clockrate_long(dat)
  dat_long$rate_diff <- TRUE
  dat_long$seed <- seed
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
  theme_bw()+
  xlab(label = "Edit rate indices")+
  ylab(label = "Posterior edit rate per day") +
  theme(legend.position = "top",
    panel.background = element_blank()
  )

p
svg(filename = paste0(figure_dir, "clock_rates.svg"), 
    width = 7.1, height = 3.1)
p
dev.off()

ggsave(filename = paste0(figure_dir, "clock_rates.png"), plot = p, dpi = 300, 
       width = 18, height = 10, units = "cm")
