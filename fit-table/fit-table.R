cat("Make a table out of fit results")

suppressPackageStartupMessages({
  library(tidyverse)
  library(kableExtra)
})

fit_dir <- here::here("fit")
fit_table_dir <- here::here("fit-table")

# Functions ===================================================================

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

make_table <- function(fits, name) {
  models <- c(
    "titre" = "Titre",
    "ili" = "ILI",
    "seroprotection" = "Seroprotection",
    "seroprotection_combined" = "Combined seroprotection",
    "seroconversion" = "seroconversion",
    "seroconversion_combined" = "Combined seroconversion"
  )
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
}

# Script ======================================================================

walk(c("titre", "seroprotection", "seroconversion"), function(name) {
  read_fits(name, exp) %>%
    select(virus, Term = term_lbl, Estimate) %>%
    pivot_wider(names_from = "virus", values_from = "Estimate") %>%
    make_table(name)
})

fits_ili <- read_fits("ili", exp) %>%
  select(Term = term_lbl, Estimate) %>%
  make_table("ili")

walk(c("seroprotection_combined", "seroconversion_combined"), function(name) {
  read_fits(name, exp) %>%
    mutate(
      n_prot = recode(
        n_prot,
        "1" = "1 antigen", "2" = "2 antigens", "3" = "3 antigens"
      )
    ) %>%
    select(n_prot, Term = term_lbl, Estimate) %>%
    pivot_wider(names_from = "n_prot", values_from = "Estimate") %>%
    make_table(name)
})
