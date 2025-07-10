library(tidyverse)
library(HDInterval)
library(ggplot2)
library(dirmult)

# set figure settings
text_size=10
label_size = 12
point_size_svg = 14
figure_dir = "sensitivity_analyses/plots/"


paths = list(
sensitivity_1 = "sensitivity_analyses/gastruloid/4-mGASv2-skyline-ou-sens1/4-mGASv2-skyline-ou-sens1.resample10000.combined.log",
sensitivity_2 = "sensitivity_analyses/gastruloid/4-mGASv2-skyline-ou-sens2/4-mGASv2-skyline-ou-sens2.resample10000.combined.log",
sensitivity_3 = "sensitivity_analyses/gastruloid/4-mGASv2-skyline-ou-sens3/4-mGASv2-skyline-ou-sens3.resample10000.combined.log",
sensitivity_4 = "sensitivity_analyses/gastruloid/4-mGASv2-skyline-ou-sens4/4-mGASv2-skyline-ou-sens4.resample10000.combined.log",
sensitivity_5 = "sensitivity_analyses/gastruloid/4-mGASv2-skyline-ou-sens5/4-mGASv2-skyline-ou-sens5.resample10000.combined.log"
)

data_list <- lapply(paths, read.table, header = TRUE)
names(data_list) <- names(paths)

# from 1 to 8
order_of_tbcs_in_alignment =c("AGGCTAATTCCC", "TAACGAAGATTT", "AACTGAATGTTT",
                              "GTATAAAGTTTG", "TTGATAACGTGA", "GTTGAAAGGTGA", 
                              "GTACAAAATAGT", "TTCACAAATTTA") 
# I want tapes 1,2,6 from Fig 5C, 
# These are the tape names from https://github.com/seidels/sciphy-materials/blob/main/figure_5/tape_seq_identifier_map.csv
# and correspond to  "AACTGAATGTTT", "AGGCTAATTCCC", "TAACGAAGATTT", which means I need the columns:
# clockRate3, clockRate1, clockRate2

# Define columns 
clock_rate_cols <- c(65, 63, 64)


process_clock_data <- function(typewriter_data, dataset_label,
                               clock_rate_cols) {
  clock_rate <- typewriter_data[, clock_rate_cols]
  
  set.seed(42)
  
  n = nrow(clock_rate)
  prior <- switch(dataset_label,
                  "sensitivity_1" = rlnorm(n, meanlog = -1.1, sdlog = 1.5),      # wider tails
                  "sensitivity_2" = rlnorm(n, meanlog = -2, sdlog = 0.8),    # shifted mean
                  "sensitivity_3" = rgamma(n, shape = 13.3, rate = 6.6),          # gamma prior
                  "sensitivity_4" = rlnorm(n, meanlog = -1.1, sdlog = 0.8),      # main text prior
                  "sensitivity_5" = rlnorm(n, meanlog = -1.1, sdlog = 0.8),      # main text prior
                  stop(paste("Unknown dataset label:", dataset_label))
  )
  clock_rate <- bind_cols(clock_rate, prior = prior)
  
  # Convert to long format and annotate source
  clock_rate_long <- clock_rate %>%
    pivot_longer(cols = everything()) %>%
    mutate(dataset = dataset_label)
  
  return(clock_rate_long)
}

processed_list <- mapply(
  process_clock_data,
  typewriter_data = data_list,
  dataset_label = names(data_list),
  MoreArgs = list(clock_rate_cols = clock_rate_cols),
  SIMPLIFY = FALSE
)
processed_list

clock_rate_all <- bind_rows(processed_list)

## only plot 95% HPD
clock_rate_hpd <- clock_rate_all %>%
  group_by(dataset, name) %>%
  summarise(
    hpd = list(hdi(value, credMass = 0.95)),
    .groups = "drop"
  ) %>%
  unnest_wider(hpd) %>%
  rename(lower = `lower`, upper = `upper`)

clock_rate_trimmed <- clock_rate_all %>%
  inner_join(clock_rate_hpd, by = c("dataset", "name"))
clock_rate_trimmed <- clock_rate_trimmed %>%
  filter(name != "prior" | value < 10)

# Plot with facet grid
colors <- c("#D95F01", "#E6AB03", "#1C9E76","#FFF2D8")
p_clock_density_facet <- ggplot(clock_rate_trimmed, aes(x = value, fill = name))+
  geom_density(alpha = 0.6, position = "identity",  bw=0.05, n=512) +
  facet_grid(rows = vars(dataset), scales = "free_y") +
  scale_fill_manual(values = colors, label = c("clockRate_3" = "Posterior Tape 1", 
                                               "clockRate_2" = "Posterior Tape 6",
                                               "clockRate_1" = "Posterior Tape 2", 
                                               "prior" = "Prior")) +
  coord_cartesian(xlim = c(0.01, 3),  expand = FALSE) + 
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 2)) +
  labs(y = "Density", x = "Edit rate ") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.80, 0.91),
    text = element_text(size = text_size),
    strip.text = element_blank(),
    axis.text.y = element_text(size = text_size -2, margin = margin(r = 0, l = -0.5)),
    axis.text.x = element_text(size = text_size - 2,
                               margin = margin(l = -0.5, r = -0.5, b = -3)),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_clock_density_facet

svg(filename = paste0(figure_dir, "gastruloid_clock_rate_sens.svg"), width = 4, height = 9,
    pointsize = point_size_svg, family = "Arial")
p_clock_density_facet
dev.off()

 ## Do stuff on the insertion probabilities
process_insert_data <- function(typewriter_data, dataset_label) {
  insert_probs <- typewriter_data[, c(15, 21, 30)]
  colnames(insert_probs) <- c("CTT", "ATA", "GAA")

  # Reorder and add prior
  ordered_names <- names(insert_probs)[order(-sapply(insert_probs, median))]
  insert_probs <- insert_probs[, ordered_names]
  set.seed(42)
  n = nrow(insert_probs)
  prior <- switch(dataset_label,
                  "sensitivity_1" = rdirichlet(nrow(insert_probs), rep(1.5,42))[,2],      # wider tails
                  "sensitivity_2" = rdirichlet(nrow(insert_probs), rep(1.5,42))[,2],    # shifted mean
                  "sensitivity_3" = rdirichlet(nrow(insert_probs), rep(1.5,42))[,2],          # gamma prior
                  "sensitivity_4" = rdirichlet(nrow(insert_probs), rep(0.5,42))[,2],      # main text prior
                  "sensitivity_5" = rdirichlet(nrow(insert_probs), rep(5,42))[,2],      # main text prior
                  stop(paste("Unknown dataset label:", dataset_label))
  )
  insert_probs <- bind_cols(insert_probs, prior = prior)
  
  # Convert to long format and annotate source
  insert_probs_long <- insert_probs %>%
    pivot_longer(cols = everything()) %>%
    mutate(name = fct_relevel(name, ordered_names, "prior"),
           dataset = dataset_label)
  
  return(insert_probs_long)
}

processed_list <- mapply(
  process_insert_data,
  typewriter_data = data_list,
  dataset_label = names(data_list),
  SIMPLIFY = FALSE
)
processed_list

insert_probs_all <- bind_rows(processed_list)

## only plot 95% HPD
insert_probs_hpd <- insert_probs_all %>%
  group_by(dataset, name) %>%
  summarise(
    hpd = list(hdi(value, credMass = 0.95)),
    .groups = "drop"
  ) %>%
  unnest_wider(hpd) %>%
  rename(lower = `lower`, upper = `upper`)

#insert_probs_trimmed <- insert_probs_all %>%
#  inner_join(insert_probs_hpd, by = c("dataset", "name")) %>%
#  filter(value >= lower, value <= upper)


# Use colours from figure 3
colors = c("#E41A1B", "#F781BF", "#000000", "#FFF2D8")
p_insert_density_facet <- ggplot(insert_probs_all, aes(x = value, fill = name))+
  geom_density(alpha = 0.6, position = "identity", bw=0.01) +
  facet_grid(rows = vars(dataset), scales = "free") +
  scale_fill_manual(values = colors, labels = c("prior" = "Prior", "CTT" = "Posterior CTT", "ATA" = "Posterior ATA", "GAA" = "Posterior GAA")) +
  coord_cartesian(xlim = c(0, 0.4), expand = FALSE) +
  #scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
  labs(y = "Density", x = "Insert probability") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.48, 0.92),
    text = element_text(size = text_size),
    strip.text = element_blank(),
    axis.text.y = element_text(size = text_size-2, margin = margin(r = 0, l = -0.5)),
    axis.text.x = element_text(size = text_size - 2,
                               margin = margin(l = -0.5, r = -0.5, b = -3)),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_insert_density_facet

svg(filename = paste0(figure_dir, "gastruloid_insert_probs_sens.svg"), width = 4, height = 9, pointsize = point_size_svg, family = "Arial")
p_insert_density_facet
dev.off()


process_growth_data <- function(typewriter_data, dataset_label) {
  growth_rates <- typewriter_data[, c("birthRate.1", "birthRate.2", "birthRate.3")] - 
    typewriter_data[, c("deathRate.1", "deathRate.2", "deathRate.3")]
  
  colnames(growth_rates) = c("growthRate.1", "growthRate.2", "growthRate.3")
  # add prior
  set.seed(42)
  n = nrow(growth_rates)
  prior = rlnorm(n = n, meanlog = -0.6, sdlog = 1) - rlnorm(n = n, meanlog = -2, sdlog = 1)
  growth_rates <- bind_cols(growth_rate = growth_rates, prior = prior)
  
  # Convert to long format 
  growth_rates <- growth_rates %>%
    pivot_longer(cols = everything()) %>%
    mutate(dataset = dataset_label)
  
  return(growth_rates)
}

processed_list <- mapply(
  process_growth_data,
  typewriter_data = data_list,
  dataset_label = names(data_list),
  SIMPLIFY = FALSE
)
processed_list
growth_rate_all = bind_rows(processed_list)

growth_rate_hpd <- growth_rate_all %>%
  group_by(dataset, name) %>%
  summarise(
    hpd = list(hdi(value, credMass = 0.95)),
    .groups = "drop"
  ) %>%
  unnest_wider(hpd) %>%
  rename(lower = `lower`, upper = `upper`)

growth_rate_trimmed <- growth_rate_all %>%
  inner_join(growth_rate_hpd, by = c("dataset", "name")) %>%
  filter(value >= lower, value <= upper)

colors = c("#d0e6f1", "#9bc2d5", "#4a90a4", "#FFF2D8")
p_growth_density <- ggplot(growth_rate_trimmed, aes(x = value, fill = name))+
  geom_density(alpha = 0.6, position = "identity", bw = 0.001) +
  facet_grid(rows = vars(dataset), scales = "free") +
  scale_fill_manual(values = colors, labels = c("prior" = "Prior",
                                                "growthRate.1" = "Posterior (Time 1)",
                                                "growthRate.2" = "Posterior (Time 2)",
                                                "growthRate.3" = "Posterior (Time 3)")) +
  coord_cartesian(xlim = c(0, 1.5), expand = FALSE) +
  #scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
  labs(y = "Density", x = "Growth rate") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.86, 0.9),
    text = element_text(size = text_size),
    strip.text = element_blank(),
    axis.text.y = element_text(size = text_size-2, margin = margin(r = 0, l = -0.5)),
    axis.text.x = element_text(size = text_size -2,
                               margin = margin(l = -0.5, r = -0.5, b = -3)),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_growth_density

svg(filename = paste0(figure_dir, "gastruloid_growth_rate_sens.svg"), width = 4, height = 9, pointsize = point_size_svg, family = "Arial")
p_growth_density
dev.off()
