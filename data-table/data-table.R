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
  group_by(virus, group, timepoint_lbl) %>%
  summarise(
    n_nomiss = sum(!is.na(titre)),
    miss = list(id[is.na(titre)]),
    .groups = "drop"
  )

# Make sure the missing are the same for all viruses
my_all_equal <- function(vec_of_4lists) {
  all(vec_of_4lists[[1]] == vec_of_4lists[[2]]) &
    all(vec_of_4lists[[2]] == vec_of_4lists[[3]]) &
    all(vec_of_4lists[[3]] == vec_of_4lists[[4]])
}
miss_counts_test <- miss_counts %>%
  group_by(group, timepoint_lbl) %>%
  summarise(all_eq = my_all_equal(miss), .groups = "drop") %>%
  filter(!all_eq)
stopifnot(nrow(miss_counts_test) == 0)

nobs <- data %>%
  filter(!is.na(titre)) %>%
  group_by(group, timepoint_lbl) %>%
  summarise(n_nomiss = length(unique(id)), .groups = "drop")

nobs %>%
  pivot_wider(names_from = "group", values_from = n_nomiss) %>%
  rename(Timepoint = timepoint_lbl) %>%
  save_csv("nobs") %>%
  kable(
    format = "latex",
    caption =
      "Number of observations at the four time points for the two groups.
      Counts apply to all four viruses.",
    label = "nobs",
    booktabs = TRUE,
    align = "lcc"
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

# Seroprotection/conversion tables

process_prot_and_conv <- function(data,
                                  group_var,
                                  count_var = "seroprotection") {
  data %>%
    filter(!is.na(!!rlang::sym(count_var))) %>%
    group_by(!!rlang::sym(group_var), group) %>%
    summarise(
      n_var = sum(!!rlang::sym(count_var)),
      n_notvar = n() - sum(!!rlang::sym(count_var)),
      prop = sum(!!rlang::sym(count_var)) / n(),
      prop_low = PropCIs::exactci(
        sum(!!rlang::sym(count_var)), n(), 0.95
      )$conf.int[[1]],
      prop_high = PropCIs::exactci(
        sum(!!rlang::sym(count_var)), n(), 0.95
      )$conf.int[[2]],
      .groups = "drop"
    ) %>%
    group_by(!!rlang::sym(group_var)) %>%
    mutate(
      `p-value` = fisher.test(
        matrix(c(n_var, n_notvar), nrow = 2)
      )$p.value
    ) %>%
    ungroup() %>%
    mutate_if(is.numeric, ~ signif(., 2)) %>%
    mutate(
      prop_est = glue::glue("{prop} ({prop_low}, {prop_high})")
    ) %>%
    select(!!rlang::sym(group_var), Group = group, prop_est, `p-value`)
}

data_seroprotection <- read_data("seroprotection")
data_seroconversion <- read_data("seroconversion")

seroprot <- process_prot_and_conv(
  data_seroprotection, "virus", "seroprotection"
)
seroconv <- process_prot_and_conv(
  data_seroconversion, "virus", "seroconversion"
)

virus_table_conv_prot <- function(data, what_prop, label, grouping = "virus") {
  data %>%
    select(-`p-value`) %>%
    pivot_wider(names_from = all_of(grouping), values_from = "prop_est") %>%
    save_csv(label) %>%
    kable(
      format = "latex",
      caption = paste(
        "Estimate (95\\% CI) of", what_prop, "proportion in the two groups.
        Confidence bounds were calculated using the Clopper-Pearson method
        as implemented in the PropCIs \\cite{PropCIs}
        R \\cite{R} package."
      ),
      label = label,
      booktabs = TRUE,
      align = "lcccc"
    ) %>%
    kable_styling(latex_options = "striped") %>%
    save_table(label)
}

virus_table_conv_prot(seroprot, "seroprotected", "prop-seroprotection")
virus_table_conv_prot(seroconv, "seroconverted", "prop-seroconversion")

virus_pvals_table_conv_prot <- function(data,
                                        what_prop,
                                        label,
                                        grouping = "virus",
                                        grouping_lab = "Virus") {
  data %>%
    pivot_wider(names_from = "Group", values_from = "prop_est") %>%
    select(
      !!rlang::sym(grouping_lab) := !!rlang::sym(grouping),
      `Standard Dose`, `High Dose`, `p-value`
    ) %>%
    save_csv(label) %>%
    kable(
      format = "latex",
      caption = paste(
        "Estimate (95\\% CI) of", what_prop, "proportion in the two groups.
        Confidence bounds were calculated using the Clopper-Pearson method
        as implemented in the PropCIs \\cite{PropCIs} R \\cite{R} package.
        The p-values were calculated using Fisher's test."
      ),
      label = label,
      booktabs = TRUE,
      align = "lcccc"
    ) %>%
    kable_styling(latex_options = "striped") %>%
    save_table(label)
}

virus_pvals_table_conv_prot(
  seroprot, "seroprotected", "prop-seroprotection-pvals"
)
virus_pvals_table_conv_prot(
  seroconv, "seroconverted", "prop-seroconversion-pvals"
)

# Combined seroprotection/seroconversion tables

data_seroprotection_combined <- read_data("seroprotection_combined")
data_seroconversion_combined <- read_data("seroconversion_combined")

recode_ag <- function(data) {
  data %>%
    mutate(
      n_prot = recode(
        n_prot,
        "1" = "1 antigen", "2" = "2 antigens", "3" = "3 antigens"
      )
    )
}

seroprot_comb <- process_prot_and_conv(
  data_seroprotection_combined, "n_prot"
) %>% recode_ag()

seroconv_comb <- process_prot_and_conv(
  data_seroconversion_combined, "n_prot"
) %>% recode_ag()

virus_table_conv_prot(
  seroprot_comb, "combined seroprotected", "prop-seroprotection_combined",
  "n_prot"
)
virus_table_conv_prot(
  seroconv_comb, "combined seroconverted", "prop-seroconversion_combined",
  "n_prot"
)

virus_pvals_table_conv_prot(
  seroprot_comb, "combined seroprotected", "prop-seroprotection_combined-pvals",
  "n_prot", "Antigens"
)
virus_pvals_table_conv_prot(
  seroconv_comb, "combined seroconverted", "prop-seroconversion_combined-pvals",
  "n_prot", "Antigens"
)

# Adverse events

adverse_events <- read_data("adverse_events")
totals <- nobs %>%
  filter(timepoint_lbl %in% c("Visit 1 (pre-vac1)", "Visit 2 (pre-vac2)")) %>%
  mutate(
    vaccine_index = str_replace(timepoint_lbl, "Visit (\\d).*", "\\1") %>%
      as.integer()
  ) %>%
  select(-timepoint_lbl)

# Counts
my_count <- function(adverse_events, variable, label) {
  adverse_events %>%
    group_by(group, vaccine_index, !!rlang::sym(variable)) %>%
    summarise(n = length(unique(id)), .groups = "drop") %>%
    mutate(count_what = label) %>%
    rename(count_category = !!rlang::sym(variable)) %>%
    inner_join(totals, by = c("group", "vaccine_index")) %>%
    mutate(n = glue::glue("{n} ({signif(n / n_nomiss * 100, 2)})"))
}

bind_rows(
  my_count(
    adverse_events %>%
      group_by(id, group, vaccine_index) %>%
      summarise(
        highest_severity = as.character(max(ordered(severity))),
        .groups = "drop"
      ),
    "highest_severity", "Highest Grade"
  ),
  my_count(adverse_events, "adverse_event", "Type")
) %>%
  pivot_wider(names_from = c("group", "vaccine_index"), values_from = "n") %>%
  select(
    count_what, count_category, `Standard Dose_1`, `High Dose_1`,
    `Standard Dose_2`, `High Dose_2`
  ) %>%
  mutate_all(~ replace_na(., "0")) %>%
  mutate(count_what = tools::toTitleCase(count_what)) %>%
  save_csv("adverse_events") %>%
  kable(
    format = "latex",
    caption =
      "Counts of unique individuals that experienced adverse events that
      fall into either a grade or a type. Percentages of corresponding groups at
      corresponding timepoints in parentheses.",
    label = "adverse_events",
    booktabs = TRUE,
    align = "lcccc",
    col.names = c(
      "", "", "Standard", "High", "Standard", "High"
    )
  ) %>%
  kable_styling(latex_options = "striped") %>%
  add_header_above(c(" " = 2, "Vaccine 1" = 2, "Vaccine 2" = 2)) %>%
  collapse_rows(1, valign = "top", latex_hline = "major") %>%
  save_table("adverse_events")

# Summarised counts
totals_summarised <- bind_rows(
  totals %>%
    mutate(vaccine_index = glue::glue("Vaccine {vaccine_index}")),
  totals %>%
    group_by(group) %>%
    summarise(n_nomiss = sum(n_nomiss), vaccine_index = "Any", .groups = "drop")
)
bind_rows(
  adverse_events %>%
    group_by(group, vaccine_index) %>%
    summarise(n = length(unique(id)), .groups = "drop") %>%
    mutate(vaccine_index = glue::glue("Vaccine {vaccine_index}")),
  adverse_events %>%
    group_by(group) %>%
    summarise(n = length(unique(id)), .groups = "drop") %>%
    mutate(vaccine_index = "Any")
) %>%
  inner_join(totals_summarised, by = c("group", "vaccine_index")) %>%
  mutate(n_noadv = n_nomiss - n) %>%
  group_by(vaccine_index) %>%
  mutate(
    `p-value` = fisher.test(matrix(c(n, n_noadv), nrow = 2))$p.value %>%
      signif(2),
    `p-value` = if_else(vaccine_index == "Any", as.character(`p-value`), ""),
    n = glue::glue("{n} ({signif(n / n_nomiss * 100, 2)})")
  ) %>%
  select(-n_nomiss, -n_noadv) %>%
  pivot_wider(names_from = "group", values_from = "n") %>%
  select(vaccine_index, `Standard Dose`, `High Dose`, `p-value`) %>%
  save_csv("adverse_events_summary") %>%
  kable(
    format = "latex",
    caption =
      "Summarised counts of unique individuals with adverse events.
      Percentages of corresponding groups at
      corresponding timepoints in parentheses. The p-value was calculated using
      Fisher's test and it relates to comparing the total number of people who
      experienced any adverse event at any time in the two groups.",
    label = "adverse_events_summary",
    booktabs = TRUE,
    align = "lccc",
    col.names = c("Timepoint", "Standard", "High", "p-value")
  ) %>%
  kable_styling(latex_options = "striped") %>%
  save_table("adverse_events_summary")
