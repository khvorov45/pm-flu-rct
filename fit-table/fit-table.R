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

exp_beta <- function(beta_name) {
  if (str_detect(beta_name, "\\{")) {
    beta_name <- str_replace(beta_name, "\\{(.*)\\}", "\\{\\\\text\\{\\1\\}\\}")
  }
  paste0("$\\text{exp}(\\beta_", beta_name, ")$")
}

save_table <- function(table_tex, table_name) {
  write(table_tex, file.path(fit_table_dir, paste0(table_name, ".tex")))
}

# Script ======================================================================

fits <- read_csv(file.path(fit_dir, "fits.csv"), col_types = cols()) %>%
  mutate(
    term_lbl = factor(
      term,
      levels = c(
        "(Intercept)",
        "groupHigh Dose",
        "timepoint_lblVisit 3", "timepoint_lblVisit 4",
        "age_years_centered", "days_since_tx_centered",
        "logtitre_baseline_centered",
        "sd_(Intercept).id", "sd_Observation.Residual"
      ),
      labels = c(
        exp_beta("0"),
        exp_beta("{HD}"),
        exp_beta("{V3}"), exp_beta("{V4}"),
        exp_beta("{AC}"), exp_beta("{XC}"),
        exp_beta("{BC}"),
        "$r_{Random}$", "$r_{Residual}$"
      )
    ),
    estimate_low = estimate - qnorm(0.975) * std.error,
    estimate_high = estimate + qnorm(0.975) * std.error,
  ) %>%
  mutate_at(vars(starts_with("estimate")), ~ cond_exp(., term_lbl)) %>%
  select(
    virus, term, term_lbl, estimate,
    std_error = std.error,
    estimate_low, estimate_high
  )

fits_ltx <- fits %>%
  select(-term) %>%
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
    Numbers in parentheses are the bounds of the 95\\% confidence interval",
    label = "estimates",
    escape = FALSE,
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(
    latex_options = "striped"
  )
save_table(fits_ltx, "fit-table")

fits_interpret <- fits %>%
  filter(virus == first(virus)) %>%
  mutate(
    Interpretation = case_when(
      term == "(Intercept)" ~ paste(
        "Expected titre for the standard dose group",
        "at visit 2, age 50, 30 days since transplant",
        "and a baseline titre measurement of 5."
      ),
      term == "groupHigh Dose" ~ paste(
        "Expected fold-titre change for the high dose group",
        "at visits 2, 3 and 4.",
        "Adjusted for age, time from transplant and baseline titre."
      ),
      term == "timepoint_lblVisit 3" ~ paste(
        "Expected fold-titre change for either group",
        "at visit 3 as compared to visit 2.",
        "Adjusted for age, time from transplant and baseline titre."
      ),
      term == "timepoint_lblVisit 4" ~ paste(
        "Expected fold-titre change for either group",
        "at visit 4 as compared to visit 2.",
        "Adjusted for age, time from transplant and baseline titre."
      ),
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
