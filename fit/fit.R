cat("Fit a model\n")

suppressPackageStartupMessages(library(tidyverse))

data_dir <- here::here("data")
fit_dir <- here::here("fit")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

fit_model <- function(data) {
  lme4::lmer(
    logtitre_mid ~ group + timepoint_lbl + age_years_centered +
      weeks4_since_tx_centered + logtitre_baseline_centered + (1 | id),
    data
  ) %>%
    broom::tidy()
}

# Script ======================================================================

data <- read_data()

data_reorg <- data %>%
  group_by(id, virus) %>%
  mutate(
    logtitre_baseline = log2(exp(logtitre[timepoint == 1L])),
    logtitre_baseline_centered = logtitre_baseline - log2(5),
  ) %>%
  ungroup() %>%
  filter(timepoint != 1L)

fits <- data_reorg %>%
  group_by(virus) %>%
  group_modify(~ fit_model(.x))
write_csv(
  fits,
  file.path(fit_dir, "fits.csv")
)
