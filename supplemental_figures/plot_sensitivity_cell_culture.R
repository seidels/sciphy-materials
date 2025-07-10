library(tidyverse)
library(HDInterval)
library(ggplot2)


# set figure settings
text_size=10
label_size = 12
point_size_svg = 14
figure_dir = "sensitivity_analyses/plots/"


paths = list(
sensitivity_1 = "sensitivity_analyses/cell_culture/1-typewriter_clockPerSite_13Sites_1000Cells_sens1/1-typewriter_clockPerSite_13Sites_1000Cells_sens1.resample10000.combined.log",
sensitivity_2 = "sensitivity_analyses/cell_culture/1-typewriter_clockPerSite_13Sites_1000Cells_sens2/1-typewriter_clockPerSite_13Sites_1000Cells_sens2.resample10000.combined.log",
sensitivity_3 = "sensitivity_analyses/cell_culture/1-typewriter_clockPerSite_13Sites_1000Cells_sens3/1-typewriter_clockPerSite_13Sites_1000Cells_sens3.resample10000.combined.log",
sensitivity_4 = "sensitivity_analyses/cell_culture/1-typewriter_clockPerSite_13Sites_1000Cells_sens4/1-typewriter_clockPerSite_13Sites_1000Cells_sens4.resample10000.combined.log",
sensitivity_5 = "sensitivity_analyses/cell_culture/1-typewriter_clockPerSite_13Sites_1000Cells_sens5/1-typewriter_clockPerSite_13Sites_1000Cells_sens5.resample10000.combined.log"
)

data_list <- lapply(paths, read.table, header = TRUE)
names(data_list) <- names(paths)

# Define shared constants
clock_rate_cols <- 28:40
col_names <- c("ATGGTAAG", "ATTTATAT", "ATTTGGTT", "GCAGGGTG", "GTAAAGAT", "TAGATTTT",
               "TGCGATTT", "TGGACGAC", "TGGTTTTG", "TTAGATTG", "TTGAGGTG", "TTTCGTGA",
               "TTCACGTA")
clock_rates_to_plot <- c("TTCACGTA", "TGGTTTTG", "TAGATTTT")

process_clock_data <- function(typewriter_data, dataset_label,
                               col_names, clock_rates_to_plot) {
  clock_rate <- typewriter_data[, 28:40]
  colnames(clock_rate) <- col_names
  clock_rate <- clock_rate %>% select(all_of(clock_rates_to_plot))
  
  # Reorder and add prior
  ordered_names <- names(clock_rate)[order(-sapply(clock_rate, median))]
  clock_rate <- clock_rate[, ordered_names]
  set.seed(42)
  
  n = nrow(clock_rate)
  prior <- switch(dataset_label,
                  "sensitivity_1" = rlnorm(n, meanlog = -2, sdlog = 1.0),      # wider tails
                  "sensitivity_2" = rlnorm(n, meanlog = -1.5, sdlog = 0.5),    # shifted mean
                  "sensitivity_3" = rgamma(n, shape = 40, rate = 20),          # gamma prior
                  "sensitivity_4" = rlnorm(n, meanlog = -2, sdlog = 0.5),      # main text prior
                  "sensitivity_5" = rlnorm(n, meanlog = -2, sdlog = 0.5),      # main text prior
                  stop(paste("Unknown dataset label:", dataset_label))
  )
  clock_rate <- bind_cols(clock_rate, prior = prior)
  
  # Convert to long format and annotate source
  clock_rate_long <- clock_rate %>%
    pivot_longer(cols = everything()) %>%
    mutate(name = fct_relevel(name, ordered_names, "prior"),
           dataset = dataset_label)
  
  return(clock_rate_long)
}

col_names <- c("ATGGTAAG", "ATTTATAT", "ATTTGGTT", "GCAGGGTG",
               "GTAAAGAT", "TAGATTTT", "TGCGATTT", "TGGACGAC",
               "TGGTTTTG", "TTAGATTG", "TTGAGGTG", "TTTCGTGA",
               "TTCACGTA")

clock_rates_to_plot <- c("TTCACGTA", "TGGTTTTG", "TAGATTTT")

processed_list <- mapply(
  process_clock_data,
  typewriter_data = data_list,
  dataset_label = names(data_list),
  MoreArgs = list(col_names = col_names, clock_rates_to_plot = clock_rates_to_plot),
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

# Plot with facet grid
colors <- c("#0A0A82", "#404085", "#A5A5C5", "#FFF2D8")
p_clock_density_facet <- ggplot(clock_rate_trimmed, aes(x = value, fill = name))+
  geom_density(alpha = 0.6, position = "identity", bw=0.01, adjust=1.5) +
  facet_grid(rows = vars(dataset), scales = "free") +
  scale_fill_manual(values = colors, labels = c("prior" = "Prior", "TTCACGTA" = "Posterior TTCACGTA",
                                                "TGGTTTTG" = "Posterior TGGTTTTG", "TAGATTTT" = "Posterior TAGATTTT")) +
  coord_cartesian(xlim = c(0.07, 2),  expand = FALSE) + 
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 2)) +
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

svg(filename = paste0(figure_dir, "clock_rate_sens.svg"), width = 4, height = 9,
    pointsize = point_size_svg, family = "Arial")
p_clock_density_facet
dev.off()

 ## Do stuff on the insertion probabilities
process_insert_data <- function(typewriter_data, dataset_label,
                               col_names, insert_probs_to_plot) {
  insert_probs <- typewriter_data[, 6:24]
  colnames(insert_probs) <- col_names
  insert_probs <- insert_probs %>% select(all_of(insert_probs_to_plot))
  
  # Reorder and add prior
  ordered_names <- names(insert_probs)[order(-sapply(insert_probs, median))]
  insert_probs <- insert_probs[, ordered_names]
  set.seed(42)
  n = nrow(insert_probs)
  prior <- switch(dataset_label,
                  "sensitivity_1" = rdirichlet(nrow(insert_probs), rep(1.5,19))[,2],      # wider tails
                  "sensitivity_2" = rdirichlet(nrow(insert_probs), rep(1.5,19))[,2],    # shifted mean
                  "sensitivity_3" = rdirichlet(nrow(insert_probs), rep(1.5,19))[,2],          # gamma prior
                  "sensitivity_4" = rdirichlet(nrow(insert_probs), rep(0.5,19))[,2],      # main text prior
                  "sensitivity_5" = rdirichlet(nrow(insert_probs), rep(5,19))[,2],      # main text prior
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

trinucleotides_names <- c("CAT","CCG","GCC","ATA","GAT","ACG","ACA","TCG","TAT","GCT","CTA","TGT","AGA","TAA","CAG","TAG","GAG","ACC","GCG")
insert_probs_to_plot = c("CAT", "GCT", "TAG")
processed_list <- mapply(
  process_insert_data,
  typewriter_data = data_list,
  dataset_label = names(data_list),
  MoreArgs = list(col_names = trinucleotides_names, insert_probs_to_plot= insert_probs_to_plot),
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

insert_probs_trimmed <- insert_probs_all %>%
  inner_join(insert_probs_hpd, by = c("dataset", "name")) %>%
  filter(value >= lower, value <= upper)

# Use colours from figure 3
colors = c("#00B814", "#7E96FF", "#E26EF7", "#FFF2D8")
p_insert_density_facet <- ggplot(insert_probs_trimmed, aes(x = value, fill = name))+
  geom_density(alpha = 0.6, position = "identity", bw=0.001) +
  facet_grid(rows = vars(dataset), scales = "free") +
  scale_fill_manual(values = colors, labels = c("CAT" = "Posterior CAT", 
                                                "prior" = "Prior", "GCT" = "Posterior GCT", "TAG" = "Posterior TAG") ) +
  coord_cartesian(xlim = c(0, 0.2), expand = FALSE) +
  #scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
  labs(y = "Density", x = "Insert probability") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.52, 0.92),
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

svg(filename = paste0(figure_dir, "insert_probs_sens.svg"), width = 4, height = 9, pointsize = point_size_svg, family = "Arial")
p_insert_density_facet
dev.off()


process_growth_data <- function(typewriter_data, dataset_label) {
  growth_rates <- typewriter_data[, "birthRate"] - typewriter_data[, "deathRate"]
  
  # add prior
  set.seed(42)
  n = length(growth_rates)
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

colors = c("white", "#FFF2D8")
p_growth_density <- ggplot(growth_rate_trimmed, aes(x = value, fill = name))+
  geom_density(alpha = 0.6, position = "identity", bw = 0.001) +
  facet_grid(rows = vars(dataset), scales = "free") +
  scale_fill_manual(values = colors, labels = c("growth_rate" = "Posterior", "prior" = "Prior")) +
  coord_cartesian(xlim = c(0, 0.75), expand = FALSE) +
  #scale_x_continuous(breaks = c(0.1, 0.2, 0.3)) +
  labs(y = "Density", x = "Growth rate") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.90, 0.9),
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

svg(filename = paste0(figure_dir, "growth_rate_sens.svg"), width = 4, height = 9, pointsize = point_size_svg, family = "Arial")
p_growth_density
dev.off()
