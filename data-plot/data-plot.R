cat("Plot the data\n")

suppressPackageStartupMessages(library(tidyverse))

data_dir <- here::here("data")
data_plot_dir <- here::here("data-plot")
data_table_dir <- here::here("data-table")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

read_table <- function(name) {
  read_csv(
    file.path(data_table_dir, paste0(name, ".csv")),
    col_types = cols()
  )
}

recode_group <- function(groups) {
  factor(
    groups,
    levels = c("Standard Dose", "High Dose"), labels = c("SD-SD", "HD-SD")
  )
}

make_plot <- function(data, data_summ) {
  data %>%
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
    scale_x_continuous("Days since first visit") +
    scale_color_discrete("Timepoint") +
    scale_shape_discrete("Timepoint") +
    facet_grid(virus ~ group) +
    geom_line(aes(group = id), alpha = 0.3) +
    geom_point(aes(shape = timepoint_lbl), alpha = 0.7, size = 1.5) +
    geom_line(
      data = data_summ,
      aes(mean_days_since_t1, titre_geom_mean), col = "#000000", lwd = 1,
      linetype = "11",
      inherit.aes = FALSE
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
}

save_plot <- function(plot, name, width, heigth) {
  ggdark::ggsave_dark(
    file.path(data_plot_dir, glue::glue("{name}.pdf")), plot,
    width = width, height = heigth, units = "cm",
  )
}

# Script ======================================================================

data <- read_data("data") %>%
  filter(!is.na(titre)) %>%
  mutate(group = recode_group(group))

data_summ <- data %>%
  group_by(virus, timepoint_lbl, group) %>%
  summarise(
    titre_geom_mean = exp(mean(logtitre_mid, na.rm = TRUE)),
    mean_days_since_t1 = mean(days_since_t1),
    .groups = "drop"
  )

# Spaghetti
spag <- make_plot(data, data_summ)
save_plot(spag, "spag", 15, 15)

spag_nobvic <- data %>%
  filter(virus != "B Vic") %>%
  make_plot(filter(data_summ, virus != "B Vic"))
save_plot(spag_nobvic, "spag-nobvic", 15, 13)

# Combined seroconversion
seroconverted_combined <- read_table("seroconversion-n_prot-long")
seroconverted_combined_plot <- seroconverted_combined %>%
  mutate(
    n_prot_lbl = ifelse(n_prot == 1, "1 strain", paste(n_prot, "strains")),
    group = recode_group(group)
  ) %>%
  ggplot(aes(group, prop)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank()
  ) +
  ylab("Proportion seroconverted") +
  ylim(c(0, 1)) +
  facet_wrap(~n_prot_lbl, strip.position = "bottom") +
  geom_errorbar(aes(ymin = prop_low, ymax = prop_high), width = 0.5) +
  geom_point(size = 4)
save_plot(seroconverted_combined_plot, "seroconverted-nprot", 12, 7)
