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
  ))

mid_est_wide <- mid_est_long %>%
  select(virus, Group = group, Timepoint = timepoint_lbl, mid_est) %>%
  mutate(Group = str_replace(Group, " Dose", "")) %>%
  pivot_wider(names_from = "virus", values_from = "mid_est")

mid_est_pval <- mid_est_long %>%
  select(group, Virus = virus, Timepoint = timepoint_lbl, mid_est) %>%
  pivot_wider(names_from = "group", values_from = "mid_est") %>%
  inner_join(
    select(
      mid_pvals,
      Virus = virus, Timepoint = timepoint_lbl, `p-value` = pval
    ),
    by = c("Virus", "Timepoint")
  )

mid_est_wide %>%
  save_csv("mid-est") %>%
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
  save_csv("mid-est-pvals") %>%
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
                                  count_var) {
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
    save_csv(glue::glue("{count_var}-{group_var}-long")) %>%
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

# Titre ratios

ratio_predata <- data_seroconversion %>%
  select(virus, id, group, contains("ratio")) %>%
  pivot_longer(contains("ratio"), names_to = "when", values_to = "ratio") %>%
  filter(!is.na(ratio))

ratio_data <- ratio_predata %>%
  group_by(virus, when, group) %>%
  summarise(
    logmean = mean(log(ratio)),
    logmean_low = logmean - qnorm(0.975) * sd(log(ratio)),
    logmean_high = logmean + qnorm(0.975) * sd(log(ratio)),
    across(contains("logmean"), ~ signif(exp(.), 2)),
    logmean_est = glue::glue("{logmean} ({logmean_low}, {logmean_high})"),
    .groups = "drop"
  )

ratio_pvals <- ratio_predata %>%
  group_by(virus, when) %>%
  summarise(
    `p-value` = t.test(
      log(ratio[group == "Standard Dose"]), log(ratio[group == "High Dose"])
    )$p.value %>% signif(3),
    .groups = "drop"
  )

ratio_data_final <- ratio_data %>%
  inner_join(ratio_pvals, by = c("virus", "when")) %>%
  mutate(
    when = str_replace(when, "ratio_", "") %>% str_replace_all("v", "V") %>%
      str_replace("_", " to ")
  ) %>%
  select(
    Virus = virus, Timepoints = when, Group = group, gmr_est = logmean_est,
    `p-value`
  )

# No p-value table
ratio_data_final %>%
  select(-`p-value`) %>%
  pivot_wider(names_from = "Virus", values_from = "gmr_est") %>%
  save_csv("gmr") %>%
  kable(
    format = "latex",
    caption = "Geometric mean titre ratios for the four viruses. 95\\% CI in
    brackets. Calculated using the normal approximation.",
    label = "gmr",
    booktabs = TRUE,
    align = "lccccc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  collapse_rows(1, valign = "top", latex_hline = "major") %>%
  save_table("gmr")

# With p-values
ratio_data_final %>%
  pivot_wider(names_from = "Group", values_from = "gmr_est") %>%
  select(Virus, Timepoints, `Standard Dose`, `High Dose`, `p-value`) %>%
  save_csv("gmr-pvals") %>%
  kable(
    format = "latex",
    caption = "Geometric mean titre ratios for the four viruses.
    95\\% CI in brackets. P-values were
    calculated with a t-test on log ratios.",
    label = "gmr-pvals",
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  collapse_rows(1, valign = "top", latex_hline = "major") %>%
  save_table("gmr-pvals")

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
  data_seroprotection_combined, "n_prot", "seroprotection"
) %>% recode_ag()

seroconv_comb <- process_prot_and_conv(
  data_seroconversion_combined, "n_prot", "seroconversion"
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
    align = "lccccc",
    col.names = c(
      "", "", "Standard group", "High group", "Standard group", "High group"
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
