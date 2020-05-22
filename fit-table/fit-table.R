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

# Script ======================================================================

fits <- read_csv(file.path(fit_dir, "fits.csv"), col_types = cols()) %>%
  mutate(
    term_lbl = factor(
      term,
      levels = c(
        "(Intercept)", "groupHigh Dose", "timepoint_lblPost-V2 Visit 1",
        "timepoint_lblPost-V2 Visit 2", "logtitre_baseline",
        "sd_(Intercept).id", "sd_Observation.Residual"
      ),
      labels = c(
        "$\\text{exp}(\\beta_0)$",
        "$\\text{exp}(\\beta_{HD})$",
        "$\\text{exp}(\\beta_{PV2-1})$",
        "$\\text{exp}(\\beta_{PV2-2})$",
        "$\\text{exp}(\\beta_{Baseline})$",
        "$r_{Random}$", "$r_{Residual}$"
      )
    ),
    estimate_low = estimate - qnorm(0.975) * std.error,
    estimate_high = estimate + qnorm(0.975) * std.error,
  ) %>%
  mutate_at(vars(starts_with("estimate")), ~ cond_exp(., term_lbl)) %>%
  select(
    virus, term_lbl, estimate,
    std_error = std.error,
    estimate_low, estimate_high
  )

fits_ltx <- fits %>%
  mutate_if(is.numeric, ~ replace_na(as.character(round(., 2)), "")) %>%
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
write(fits_ltx, file.path(fit_table_dir, "fit-table.tex"))
