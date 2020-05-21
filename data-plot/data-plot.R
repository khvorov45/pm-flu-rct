# Plots of the data

library(tidyverse)

data_dir <- here::here("data")
data_plot_dir <- here::here("data-plot")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

# Script ======================================================================

data <- read_data() %>% filter(!is.na(titre))

data_summ <- data %>%
  group_by(virus, timepoint_lbl, group) %>%
  summarise(
    titre_geom_mean = exp(mean(log(titre), na.rm = TRUE)),
    mean_days_since_t1 = mean(days_since_t1)
  )

spag <- data %>%
  filter(!is.na(titre)) %>%
  ggplot(aes(days_since_t1, titre, col = timepoint_lbl)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
  ) +
  scale_y_log10("HI Titre", breaks = 5 * 2^(0:11)) +
  scale_x_continuous("Days since T1") +
  scale_color_discrete("Timepoint") +
  facet_grid(group ~ virus) +
  geom_line(aes(group = id), alpha = 0.3) +
  geom_point(aes(col = timepoint_lbl)) +
  geom_line(
    data = data_summ,
    aes(mean_days_since_t1, titre_geom_mean), col = "#000000", lwd = 1,
    linetype = "11",
    inherit.aes = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggdark::ggsave_dark(
  file.path(data_plot_dir, "spag.pdf"), spag,
  width = 20, height = 10, units = "cm",
)
