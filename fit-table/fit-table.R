cat("Make a table out of fit results")

suppressPackageStartupMessages({
  library(tidyverse)
  library(kableExtra)
})

fit_dir <- here::here("fit")
fit_table_dir <- here::here("fit-table")

# Functions ===================================================================

cond_exp <- function(estimate, term) {
  estimate <- if_else(str_detect(term, "^\\$r_"), estimate, exp(estimate))
}

gen_vars_ltx <- function(varsint) {
  pwalk()
}

save_table <- function(table_tex, table_name) {
  write(table_tex, file.path(fit_table_dir, paste0(table_name, ".tex")))
}

# Script ======================================================================

fits <- read_csv(file.path(fit_dir, "fits.csv"), col_types = cols()) %>%
  mutate(
    estimate_low = estimate - qnorm(0.975) * std.error,
    estimate_high = estimate + qnorm(0.975) * std.error,
  ) %>%
  mutate_at(vars(starts_with("estimate")), ~ cond_exp(., term_lbl)) %>%
  select(
    virus, term, term_lbl, var_lbl, estimate,
    std_error = std.error,
    estimate_low, estimate_high
  )

fits_ltx <- fits %>%
  mutate_if(is.numeric, ~ replace_na(as.character(signif(., 2)), "")) %>%
  mutate(
    Estimate = glue::glue("{estimate} ({estimate_low}, {estimate_high})") %>%
      str_replace(" \\(, \\)", "")
  ) %>%
  select(virus, Term = term_lbl, Estimate) %>%
  pivot_wider(names_from = "virus", values_from = "Estimate") %>%
  kable(
    format = "latex",
    caption = "Model parameter estimates for the four viruses.
    Numbers in parentheses are the bounds of the 95\\% confidence interval.",
    label = "estimates",
    escape = FALSE,
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(
    latex_options = "striped"
  )
save_table(fits_ltx, "fit-table")

vars_tbl <- fits %>%
  filter(!is.na(var_lbl)) %>%
  select(var_lbl) %>%
  distinct() %>%
  mutate(
    interpret = case_when(
      var_lbl == "$T_{l,m}$" ~ "Natural log-titre shifted to the midpoint",
      var_lbl == "$G_H$" ~ paste(
        "Indicator of whether the observation belongs to the",
        "high-dose group (1) or the standard-dose group (0)."
      ),
      var_lbl == "$V_3$" ~ paste(
        "Indicator of whether the observation belongs to the",
        "visit 3 (1) or not (0)."
      ),
      var_lbl == "$V_4$" ~ paste(
        "Indicator of whether the observation belongs to the",
        "visit 4 (1) or not (0)."
      ),
      var_lbl == "$M$" ~
      "Indicator of whether the subject has myeloma (1) or not (0)",
      var_lbl == "$P$" ~ paste(
        "Indicator of whether the subject was vaccinated no earlier than 2018",
        "(1) or either vaccinated earlier than 2018 or never vaccinated (0)."
      ),
      var_lbl == "$A_C$" ~ "Age in years at first visit. Centred on 50.",
      var_lbl == "$X_C$" ~ paste(
        "Number of 4-week periods (approximately months) from transplant",
        "to first visit. Centred on 1."
      ),
      var_lbl == "$B_C$" ~ paste(
        "Log2-baseline titre (i.e. titre at visit 1) shifted to midpoint.",
        "Centred on $\\text{log}_2(5)$."
      ),
      TRUE ~ "Some other variable"
    ),
    varline = paste(
      var_lbl, "---", interpret, ifelse(
        row_number() == max(row_number()), "", "\\\\"
      )
    )
  )
writeLines(vars_tbl$varline, file.path(fit_table_dir, "vars.tex"))

fits_interpret <- fits %>%
  filter(virus == first(virus)) %>%
  mutate(
    Interpretation = case_when(
      term == "(Intercept)" ~ paste(
        "Expected titre for the standard dose group",
        "at visit 2, age 50, 4 weeks since transplant,",
        "baseline titre measurement of 5 and cancer type other than myeloma."
      ),
      term == "groupHigh Dose" ~ paste(
        "Expected fold-titre change for the high dose group",
        "as compared to the standard dose group",
        "at visits 2, 3 and 4.",
        "Adjusted for age, time from transplant,",
        "baseline titre and myeloma status"
      ),
      term == "timepoint_lblVisit 3" ~ paste(
        "Expected fold-titre change for either group",
        "at visit 3 as compared to visit 2.",
        "Adjusted for age, time from transplant,",
        "baseline titre and myeloma status"
      ),
      term == "timepoint_lblVisit 4" ~ paste(
        "Expected fold-titre change for either group",
        "at visit 4 as compared to visit 2.",
        "Adjusted for age, time from transplant,",
        "baseline titre and myeloma status"
      ),
      term %in% c("age_years_centered", "age_years_baseline_centered") ~ paste(
        "Expected fold-titre increase for either group at visits 2, 3 and 4",
        "for 1 year increase in age.",
        "Adjusted for time from transplant, baseline titre and myeloma status."
      ),
      term %in% c(
        "weeks4_since_tx_centered", "weeks4_since_tx_baseline_centered"
      ) ~ paste(
        "Expected fold-titre increase for either group at visits 2, 3 and 4",
        "for a 4-week increase in time from transplant.",
        "Adjusted for age, baseline titre and myeloma status."
      ),
      term == "logtitre_baseline_centered" ~ paste(
        "Expected fold-titre increase for either group at visits 2, 3 and 4",
        "for a 2-fold increase in the baseline titre.",
        "Adjusted for age, time from transplant and myeloma status."
      ),
      term == "myeloma" ~ paste(
        "Expected fold-titre increase for either group at visits 2, 3 and 4",
        "for subjects with myeloma compared to subjects with other cancer type",
        "Adjusted for age and time from transplant."
      ),
      term == "vac_in_prior_year" ~ paste(
        "Expected fold-titre increase for either group at visits 2, 3 and 4",
        "for subjects last vaccinated no earlier than 2018",
        "compared to subjects last vaccinated earlier",
        "than 2018 or never vaccinated.",
        "Adjusted for age and time from transplant."
      ),
      term == "sd_(Intercept).id" ~ "A measure of between-subject variability.",
      term == "sd_Observation.Residual" ~
      "A measure of within-subject variability.",
      TRUE ~ "Some other interpretation"
    ),
  ) %>%
  select(Term = term_lbl, Interpretation) %>%
  kable(
    format = "latex",
    caption = "Interpretation of the model parameters",
    label = "estimates-interpret",
    escape = FALSE,
    booktabs = TRUE,
    align = "ll"
  ) %>%
  kable_styling(
    latex_options = "striped"
  ) %>%
  column_spec(2, width = "35em")
save_table(fits_interpret, "fit-interpret")
