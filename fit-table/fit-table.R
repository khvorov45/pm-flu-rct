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

fits %>%
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
  ) %>%
  save_table("fit-table")
