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

data <- read_data("data")

miss_counts <- data %>%
  group_by(virus, timepoint_lbl) %>%
  summarise(
    n_nomiss = sum(!is.na(titre)),
    .groups = "drop"
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

mid_pvals <- data %>%
  filter(!is.na(titre)) %>%
  group_by(virus, timepoint_lbl) %>%
  summarise(
    pval = t.test(
      logtitre_mid[group == "Standard Dose"],
      logtitre_mid[group == "High Dose"],
      var.equal = TRUE,
    )$p.value,
    .groups = "drop"
  )

mid_est_long <- data %>%
  filter(!is.na(titre)) %>%
  group_by(group, virus, timepoint_lbl) %>%
  summarise(
    mid_mean = mean(logtitre_mid),
    logmid_sd = sd(logtitre_mid),
    nobs = n(),
    mid_mean_lb = mid_mean - qnorm(0.975) * logmid_sd / sqrt(nobs),
    mid_mean_ub = mid_mean + qnorm(0.975) * logmid_sd / sqrt(nobs),
    .groups = "drop"
  ) %>%
  mutate(mid_est = glue::glue(
    "{signif(exp(mid_mean), 2)} ",
    "({signif(exp(mid_mean_lb), 2)}, {signif(exp(mid_mean_ub), 2)})"
  )) %>%
  save_csv("mid-long")

mid_est_wide <- mid_est_long %>%
  select(virus, Group = group, Timepoint = timepoint_lbl, mid_est) %>%
  mutate(Group = str_replace(Group, " Dose", "")) %>%
  pivot_wider(names_from = "virus", values_from = "mid_est") %>%
  save_csv("mid-wide")

mid_est_pval <- mid_est_long %>%
  select(group, Virus = virus, Timepoint = timepoint_lbl, mid_est) %>%
  pivot_wider(names_from = "group", values_from = "mid_est") %>%
  inner_join(
    select(
      mid_pvals,
      Virus = virus, Timepoint = timepoint_lbl, `p-value` = pval
    ),
    by = c("Virus", "Timepoint")
  ) %>%
  save_csv("mid-pvals")

mid_est_wide %>%
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

mid_est_pval %>%
  mutate(`p-value` = signif(`p-value`, 2)) %>%
  kable(
    format = "latex",
    caption =
      "Estimate (95\\% CI) of the geometric mean of the HI titres at the four
      timepoints for the four viruses for the two groups.
      The mean (and interval) were calculated using titre midpoints on the
      log-scale and then exponentiated. The p-value was calculated using the
      two-sample t-test with pooled variance.",
    label = "mid-est-pvals",
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  collapse_rows(columns = 1, valign = "top", latex_hline = "major") %>%
  save_table("mid-est-pvals")

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
    .groups = "drop"
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

# Seroprotection table

data_seroprotection <- read_data("seroprotection")

data_seroprotection %>%
  filter(!is.na(seroprotection)) %>%
  group_by(virus, group) %>%
  summarise(
    prop = sum(seroprotection) / n(),
    prop_low = PropCIs::exactci(sum(seroprotection), n(), 0.95)$conf.int[[1]],
    prop_high = PropCIs::exactci(sum(seroprotection), n(), 0.95)$conf.int[[2]],
    .groups = "drop"
  ) %>%
  mutate_if(is.numeric, ~ signif(., 2)) %>%
  mutate(
    prop_est = glue::glue("{prop} ({prop_low}, {prop_high})")
  ) %>%
  select(virus, Group = group, prop_est) %>%
  pivot_wider(names_from = "virus", values_from = "prop_est") %>%
  save_csv("prop-seroprotection") %>%
  kable(
    format = "latex",
    caption =
      "Estimate (95\\% CI) of seroprotected proportion in the two groups.
      Confidence bounds were calculated using the Clopper-Pearson method
      as implemented in the PropCIs \\cite{PropCIs} R \\cite{R} package.",
    label = "prop-seroprotection",
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  save_table("prop-seroprotection")

# Combined seroprotection table

data_seroprotection_combined <- read_data("seroprotection_combined")

data_seroprotection_combined %>%
  filter(!is.na(seroprotection)) %>%
  group_by(n_prot, group) %>%
  summarise(
    prop = sum(seroprotection) / n(),
    prop_low = PropCIs::exactci(sum(seroprotection), n(), 0.95)$conf.int[[1]],
    prop_high = PropCIs::exactci(sum(seroprotection), n(), 0.95)$conf.int[[2]],
    .groups = "drop"
  ) %>%
  mutate_if(is.numeric, ~ signif(., 2)) %>%
  mutate(
    prop_est = glue::glue("{prop} ({prop_low}, {prop_high})")
  ) %>%
  select(n_prot, Group = group, prop_est) %>%
  mutate(
    n_prot = recode(
      n_prot,
      "1" = "1 antigen", "2" = "2 antigens", "3" = "3 antigens"
    )
  ) %>%
  pivot_wider(names_from = "n_prot", values_from = "prop_est") %>%
  save_csv("prop-seroprotection_combined") %>%
  kable(
    format = "latex",
    caption =
      "Estimate (95\\% CI) of combined seroprotected
      proportion in the two groups.
      Confidence bounds were calculated using the Clopper-Pearson method
      as implemented in the PropCIs \\cite{PropCIs} R \\cite{R} package.",
    label = "prop-seroprotection_combined",
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  save_table("prop-seroprotection_combined")
