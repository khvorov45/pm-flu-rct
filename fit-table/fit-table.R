cat("Make a table out of fit results")

suppressPackageStartupMessages({
  library(tidyverse)
  library(kableExtra)
})

data_dir <- here::here("data")
fit_dir <- here::here("fit")
fit_table_dir <- here::here("fit-table")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

invlogit <- function(x) 1 - 1 / (1 + exp(x))

read_fits <- function(name, fun) {
  read_csv(
    file.path(fit_dir, glue::glue("fits-{name}.csv")),
    col_types = cols()
  ) %>%
    mutate(
      estimate_low = estimate - qnorm(0.975) * std.error,
      estimate_high = estimate + qnorm(0.975) * std.error,
    ) %>%
    mutate_at(vars(starts_with("estimate")), ~ cond_exp(., term_lbl, fun)) %>%
    mutate_if(is.numeric, ~ replace_na(as.character(signif(., 2)), "")) %>%
    mutate(
      Estimate = glue::glue("{estimate} ({estimate_low}, {estimate_high})") %>%
        str_replace(" \\(, \\)", "")
    )
}

cond_exp <- function(estimate, term, fun) {
  estimate <- if_else(str_detect(term, "^\\$r_"), estimate, fun(estimate))
}

save_csv <- function(data, name) {
  write_csv(data, file.path(fit_table_dir, paste0(name, ".csv")))
}

save_table <- function(table_tex, table_name) {
  write(table_tex, file.path(fit_table_dir, paste0(table_name, ".tex")))
}

make_table <- function(fits, name, models) {
  fits %>%
    save_csv(glue::glue("fit-table-{name}")) %>%
    kable(
      format = "latex",
      caption = glue::glue(
        "{models[[name]]} model parameter estimates. ",
        "Numbers in parentheses are the bounds of the 95\\% confidence interval."
      ),
      label = glue::glue("estimates-{name}"),
      escape = FALSE,
      booktabs = TRUE,
      align = "lcccc"
    ) %>%
    kable_styling(
      latex_options = "striped"
    ) %>%
    save_table(glue::glue("fit-table-{name}"))
  fits
}

make_table_pval <- function(fits, name, models) {
  tabl <- fits %>%
    save_csv(glue::glue("fit-table-{name}-pval")) %>%
    kable(
      format = "latex",
      caption = glue::glue(
        "{models[[name]]} model parameter estimates. ",
        "Numbers in parentheses are the bounds of the 95\\% confidence interval."
      ),
      label = glue::glue("estimates-{name}-pval"),
      escape = FALSE,
      booktabs = TRUE,
      align = "llcc"
    ) %>%
    kable_styling(
      latex_options = "striped"
    )
  if (names(fits)[[1]] != "Term") {
    tabl <- collapse_rows(tabl, 1, valign = "top", latex_hline = "major")
  }
  save_table(tabl, glue::glue("fit-table-{name}-pval"))
  fits
}

# Script ======================================================================

models <- c(
  "titre" = "Titre",
  "ili" = "ILI",
  "seroprotection" = "Seroprotection",
  "seroprotection_combined" = "Combined seroprotection",
  "seroconversion" = "Seroconversion",
  "seroconversion_combined" = "Combined seroconversion"
)

walk(c("titre", "seroprotection", "seroconversion"), function(name) {
  read_fits(name, exp) %>%
    select(Virus = virus, Term = term_lbl, Estimate, `p-value` = p.value) %>%
    make_table_pval(name, models) %>%
    select(-`p-value`) %>%
    pivot_wider(names_from = "Virus", values_from = "Estimate") %>%
    make_table(name, models)
})

fits_ili <- read_fits("ili", exp) %>%
  select(Term = term_lbl, Estimate, `p-value` = p.value) %>%
  make_table_pval("ili", models) %>%
  select(-`p-value`) %>%
  make_table("ili", models)

walk(c("seroprotection_combined", "seroconversion_combined"), function(name) {
  read_fits(name, exp) %>%
    mutate(
      n_prot = recode(
        n_prot,
        "1" = "1 antigen", "2" = "2 antigens", "3" = "3 antigens"
      )
    ) %>%
    select(Antigens = n_prot, Term = term_lbl, Estimate, `p-value` = p.value) %>%
    make_table_pval(name, models) %>%
    select(-`p-value`) %>%
    pivot_wider(names_from = "Antigens", values_from = "Estimate") %>%
    make_table(name, models)
})

fits_seroconversion_combined3 <- read_fits(
  "seroconversion_combined", exp
) %>%
  filter(n_prot == "3") %>%
  select(term, Estimate_mult = Estimate, `p-value_mult` = p.value)

fits_seroconversion_combined3_ind_mult <- read_fits(
  "seroconversion_combined3_individual", exp
) %>%
  mutate(
    variable_lbl = recode(
      variable_name,
      "group" = "High Dose",
      "vac_in_prior_year" = "Vaccinated in prior year",
      "current_therapy" = "Receives therapy",
      "age_years_baseline_centered" = "Age, years",
      "weeks4_since_tx_baseline_centered" =
        "Weeks from transplant"
    )
  ) %>%
  select(
    term,
    Term = variable_lbl,
    Estimate_ind = Estimate, `p-value_ind` = p.value,
    contains("eroconvert")
  ) %>%
  inner_join(
    fits_seroconversion_combined3,
    by = "term"
  ) %>%
  select(
    Term, contains("Seroconverted"), contains("Did not seroconvert"),
    Estimate_ind, `p-value_ind`, Estimate_mult, `p-value_mult`,
  )

fits_seroconversion_combined3_ind_mult %>%
  save_csv(glue::glue("fit-table-seroconversion_combined3_ind_mult")) %>%
  kable(
    format = "latex",
    caption = glue::glue(
      "Data summaries (count (proportion) or mean $\\pm$ sd),
      parameter estimates from univariate and multivariate models for
      seroconversion against 3 antigens."
    ),
    label = glue::glue("estimates-seroconversion_combined3_ind_mult"),
    escape = FALSE,
    booktabs = TRUE,
    align = "l",
    col.names = c(
      names(fits_seroconversion_combined3_ind_mult)[1:3],
      "Estimate", "p-value", "Estimate", "p-value"
    )
  ) %>%
  kable_styling(
    latex_options = c("striped", "scale_down")
  ) %>%
  add_header_above(c(" " = 3, "Univariate" = 2, "Multivariate" = 2)) %>%
  column_spec(1, width = "7em") %>%
  column_spec(2, width = "7em") %>%
  column_spec(3, width = "7em") %>%
  save_table(glue::glue("fit-table-seroconversion_combined3_ind_mult"))
