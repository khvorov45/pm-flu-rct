# Plots of the data

library(tidyverse)

data_dir <- here::here("data")
data_plot_dir <- here::here("data-plot")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

# Script ======================================================================

data <- read_data()

data_summ <- data %>%
  group_by(virus, timepoint, group) %>%
  summarise(titre_geom_mean = exp(mean(log(titre), na.rm = TRUE)))

spag <- data %>%
  filter(!is.na(titre)) %>%
  ggplot(aes(timepoint, titre)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_log10("HI Titre", breaks = 5 * 2^(0:11)) +
  scale_x_continuous("Timepoint") +
  facet_grid(group ~ virus) +
  geom_line(aes(group = id), alpha = 0.3) +
  geom_line(
    data = data_summ,
    aes(y = titre_geom_mean), col = "red", lwd = 1
  )
ggdark::ggsave_dark(
  file.path(data_plot_dir, "spag.pdf"), spag,
  width = 20, height = 10, units = "cm",
)
