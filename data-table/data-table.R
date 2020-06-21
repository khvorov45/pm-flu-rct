cat("Make data table")

suppressPackageStartupMessages({
  library(tidyverse)
  library(kableExtra)
})

data_dir <- here::here("data")
data_table_dir <- here::here("data-table")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

save_csv <- function(table, name) {
  write_csv(table, file.path(data_table_dir, paste0(name, ".csv")))
  table
}

save_table <- function(table, name) {
  write(table, file.path(data_table_dir, paste0(name, ".tex")))
}

# Script ======================================================================

data <- read_data()

miss_counts <- data %>%
  group_by(virus, timepoint_lbl) %>%
  summarise(
    n_nomiss = sum(!is.na(titre))
  )

miss_counts_tbl <- miss_counts %>%
  pivot_wider(names_from = "virus", values_from = n_nomiss) %>%
  rename(Timepoint = timepoint_lbl) %>%
  save_csv("nobs") %>%
  kable(
    format = "latex",
    caption =
      "Number of observations at the four time points for the four viruses.",
    label = "nobs",
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  save_table("nobs")

mid_est_long <- data %>%
  filter(!is.na(titre)) %>%
  group_by(group, virus, timepoint_lbl) %>%
  summarise(
    mid_mean = mean(logtitre_mid),
    logmid_sd = sd(logtitre_mid),
    nobs = n(),
    mid_mean_lb = mid_mean - qnorm(0.975) * logmid_sd / sqrt(nobs),
    mid_mean_ub = mid_mean + qnorm(0.975) * logmid_sd / sqrt(nobs),
  ) %>%
  ungroup() %>%
  save_csv("mid-long")

mid_est_long %>%
  mutate_if(is.numeric, ~ signif(exp(.), 2)) %>%
  mutate(
    mid_est = glue::glue("{mid_mean} ({mid_mean_lb}, {mid_mean_ub})")
  ) %>%
  select(virus, Group = group, Timepoint = timepoint_lbl, mid_est) %>%
  mutate(Group = str_replace(Group, " Dose", "")) %>%
  pivot_wider(names_from = "virus", values_from = "mid_est") %>%
  save_csv("mid-wide")

mid_est %>%
  kable(
    format = "latex",
    caption =
      "Estimate (95\\% CI) of the geometric mean of the HI titres at the four
      timepoints for the four viruses for the two groups.
      The mean (and interval) were calculated using titre midpoints on the
      log-scale and then exponentiated.",
    label = "mid-est",
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  collapse_rows(columns = 1, valign = "top", latex_hline = "major") %>%
  save_table("mid-est")

data_wide <- data %>%
  filter(timepoint == 1L) %>%
  pivot_wider(names_from = "virus", values_from = contains("titre"))

# Make sure there is one row per individual
stopifnot(all(data_wide$id == unique(data$id)))

# ILI proportion table

data_wide %>%
  filter(!is.na(ili)) %>%
  group_by(group) %>%
  summarise(
    prop_ili = sum(ili) / n(),
    se_prop = sqrt(prop_ili * (1 - prop_ili) / n()),
    prop_ili_low = PropCIs::exactci(sum(ili), n(), 0.95)$conf.int[[1]],
    prop_ili_high = PropCIs::exactci(sum(ili), n(), 0.95)$conf.int[[2]],
  ) %>%
  mutate_if(is.numeric, ~ signif(., 2)) %>%
  mutate(
    prop_ili_est = glue::glue("{prop_ili} ({prop_ili_low}, {prop_ili_high})")
  ) %>%
  select(Group = group, `ILI proportion` = prop_ili_est) %>%
  save_csv("prop-ili") %>%
  kable(
    format = "latex",
    caption =
      "Estimate (95\\% CI) of ILI proportion in the two groups.
      Confidence bounds were calculated using the Clopper-Pearson method
      as implemented in the PropCIs \\cite{PropCIs} R \\cite{R} package.",
    label = "prop-ili",
    booktabs = TRUE,
    align = "lc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  save_table("prop-ili")
