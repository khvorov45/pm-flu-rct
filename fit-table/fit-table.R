cat("Make a table out of fit results")

suppressPackageStartupMessages({
  library(tidyverse)
  library(kableExtra)
})

fit_dir <- here::here("fit")
fit_table_dir <- here::here("fit-table")

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
        "$\\beta_0$", "$\\beta_{HD}$", "$\\beta_{PV2-1}$",
        "$\\beta_{PV2-2}$", "$\\beta_{Baseline}$",
        "$r_{Random}$", "$r_{Residual}$"
      )
    ),
    estimate_low = estimate - qnorm(0.975) * std.error,
    estimate_high = estimate + qnorm(0.975) * std.error,
  ) %>%
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
    caption = "Model parameter estimates for the four viruses",
    label = "estimates",
    escape = FALSE
  ) %>%
  kable_styling()
write(fits_ltx, file.path(fit_table_dir, "fit-table.tex"))
