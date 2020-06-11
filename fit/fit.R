cat("Fit a model\n")

suppressPackageStartupMessages(library(tidyverse))

data_dir <- here::here("data")
fit_dir <- here::here("fit")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

exp_beta <- function(beta_name) {
  if (str_detect(beta_name, "\\{")) {
    beta_name <- str_replace(beta_name, "\\{(.*)\\}", "\\{\\\\text\\{\\1\\}\\}")
  }
  paste0("$\\text{exp}(\\beta_", beta_name, ")$")
}

fit_model <- function(data) {
  lme4::lmer(
    logtitre_mid ~ group + timepoint_lbl + myeloma +
      vac_in_prior_year +
      age_years_baseline_centered +
      weeks4_since_tx_baseline_centered + logtitre_baseline_centered + (1 | id),
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
    age_years_baseline_centered = age_years_centered[timepoint == 1L],
    weeks4_since_tx_baseline_centered =
      weeks4_since_tx_centered[timepoint == 1L],
  ) %>%
  ungroup() %>%
  filter(timepoint != 1L)

fits <- data_reorg %>%
  group_by(virus) %>%
  group_modify(~ fit_model(.x))

fits_ref <- tribble(
  ~term, ~term_lbl, ~var_lbl,
  "(Intercept)", exp_beta("0"), "$T_{l,m}$",
  "groupHigh Dose", exp_beta("{HD}"), "$G_H$",
  "timepoint_lblVisit 3", exp_beta("{V3}"), "$V_3$",
  "timepoint_lblVisit 4", exp_beta("{V4}"), "$V_4$",
  "myeloma", exp_beta("M"), "$M$",
  "vac_in_prior_year", exp_beta("{PV}"), "$P$",
  "age_years_centered", exp_beta("{AC}"), "$A_C$",
  "age_years_baseline_centered", exp_beta("{AC}"), "$A_C$",
  "weeks4_since_tx_centered", exp_beta("{XC}"), "$X_C$",
  "weeks4_since_tx_baseline_centered", exp_beta("{XC}"), "$X_C$",
  "logtitre_baseline_centered", exp_beta("{BC}"), "$B_C$",
  "sd_(Intercept).id", "$r_{Random}$", "",
  "sd_Observation.Residual", "$r_{Residual}$", "",
)

write_csv(
  left_join(fits, fits_ref, by = "term"),
  file.path(fit_dir, "fits.csv")
)
