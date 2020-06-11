cat("Fit a model\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(kableExtra)
})

data_dir <- here::here("data")
fit_dir <- here::here("fit")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

exp_beta <- function(beta_name) {
  if (str_detect(beta_name, "\\{")) {
    beta_name <- str_replace(beta_name, "\\{(.*)\\}", "\\{\\\\text\\{\\1\\}\\}")
  }
  paste0("$\\text{exp}(\\beta_", beta_name, ")$")
}

fit_model <- function(data) {
  lme4::lmer(
    logtitre_mid ~ group +
      timepoint_lbl
      + myeloma
      + vac_in_prior_year
      + current_therapy
      + age_years_baseline_centered
      + weeks4_since_tx_baseline_centered
      + logtitre_baseline_centered
      + (1 | id),
    data
  ) %>%
    broom::tidy()
}

gen_b0_int <- function() {
  all <- list(
    "myeloma" = "cancer other than myeloma",
    "vac_in_prior_year" = "never vaccinated or vaccinated earlier than 2018",
    "current_therapy" = "not receiving therapy",
    "age_years_centered" = "age 50",
    "age_years_baseline_centered" = "age 50",
    "weeks4_since_tx_centered" = "4 weeks since transplant",
    "weeks4_since_tx_baseline_centered" = "4 weeks since transplant",
    "logtitre_baseline_centered" = "baseline titre of 5"
  )
  all <- all[names(all) %in% fits$term]
  paste0(
    "Expected titre for the standard dose group at visit 2, ",
    paste(all, collapse = ", "), "."
  )
}

adjusted <- function(this = NULL) {
  adj <- list(
    "myeloma" = "myeloma status",
    "vac_in_prior_year" = "previous vaccination status",
    "current_therapy" = "current therapy status",
    "age_years_centered" = "age",
    "age_years_baseline_centered" = "age",
    "weeks4_since_tx_centered" = "time from transplant",
    "weeks4_since_tx_baseline_centered" = "time from transplant",
    "logtitre_baseline_centered" = "baseline titre"
  )
  adj <- adj[names(adj) %in% fits$term]
  if (!is.null(this)) adj <- adj[!names(adj) %in% this]
  paste0("Adjusted for ", paste(adj, collapse = ", "), ".")
}

clean_label <- function(labels) {
  str_replace_all(labels, "\\$", "") %>%
    str_replace_all("\\\\text\\{exp\\}\\((.*)\\)", "\\1")
}

# Script ======================================================================

data <- read_data()

data_reorg <- data %>%
  group_by(id, virus) %>%
  mutate(
    logtitre_baseline = log2(exp(logtitre[timepoint == 1L])),
    logtitre_baseline_centered = logtitre_baseline - log2(5),
    age_years_baseline_centered = age_years_centered[timepoint == 1L],
    weeks4_since_tx_baseline_centered =
      weeks4_since_tx_centered[timepoint == 1L],
  ) %>%
  ungroup() %>%
  filter(timepoint != 1L)

fits <- data_reorg %>%
  group_by(virus) %>%
  group_modify(~ fit_model(.x))

fits_ref <- tribble(
  ~term, ~term_lbl, ~var_lbl, ~term_int, ~var_int,

  "(Intercept)", exp_beta("0"), "$T_{l,m}$",
  gen_b0_int(), "Natural log-titre shifted to the midpoint.",

  "groupHigh Dose", exp_beta("{HD}"), "$G_H$",
  paste(
    "Expected fold-titre change for the high dose group",
    "as compared to the standard dose group",
    "at visits 2, 3 and 4.",
    adjusted()
  ),
  paste(
    "Indicator of whether the observation belongs to the",
    "high-dose group (1) or the standard-dose group (0)."
  ),

  "timepoint_lblVisit 3", exp_beta("{V3}"), "$V_3$",
  paste(
    "Expected fold-titre change for either group",
    "at visit 3 as compared to visit 2.",
    adjusted()
  ),
  paste(
    "Indicator of whether the observation belongs to",
    "visit 3 (1) or not (0)."
  ),

  "timepoint_lblVisit 4", exp_beta("{V4}"), "$V_4$",
  paste(
    "Expected fold-titre change for either group",
    "at visit 4 as compared to visit 2.",
    adjusted()
  ),
  paste(
    "Indicator of whether the observation belongs to",
    "visit 4 (1) or not (0)."
  ),

  "myeloma", exp_beta("M"), "$M$",
  paste(
    "Expected fold-titre increase for either group at visits 2, 3 and 4",
    "for subjects with myeloma compared to subjects with other cancer type",
    adjusted("myeloma")
  ),
  "Indicator of whether the subject has myeloma (1) or not (0)",

  "vac_in_prior_year", exp_beta("{PV}"), "$P$",
  paste(
    "Expected fold-titre increase for either group at visits 2, 3 and 4",
    "for subjects last vaccinated no earlier than 2018",
    "compared to subjects last vaccinated earlier",
    "than 2018 or never vaccinated.",
    adjusted("vac_in_prior_year")
  ),
  paste(
    "Indicator of whether the subject was vaccinated no earlier than 2018",
    "(1) or either vaccinated earlier than 2018 or never vaccinated (0)."
  ),

  "current_therapy", exp_beta("{R}"), "$R$",
  paste(
    "Expected fold-titre increase for either group at visits 2, 3 and 4",
    "for subjects receiving therapy as compared to subjects not",
    "receiving therapy.",
    adjusted("current_therapy")
  ),
  paste(
    "Indicator of whether the subject was receiving therapy (1)",
    "or not (0)."
  ),

  "age_years_centered", exp_beta("{AC}"), "$A_C$",
  paste(
    "Expected fold-titre increase for either group at visits 2, 3 and 4",
    "for 1 year increase in age.",
    adjusted("age_years_centered")
  ),
  "Age in years at first visit. Centred on 50.",

  "age_years_baseline_centered", exp_beta("{AC}"), "$A_C$",
  paste(
    "Expected fold-titre increase for either group at visits 2, 3 and 4",
    "for 1 year increase in age.",
    adjusted("age_years_baseline_centered")
  ),
  "Age in years at time of measurement. Centred on 50.",

  "weeks4_since_tx_centered", exp_beta("{XC}"), "$X_C$",
  paste(
    "Expected fold-titre increase for either group at visits 2, 3 and 4",
    "for a 4-week increase in time from transplant.",
    adjusted("weeks4_since_tx_centered")
  ),
  paste(
    "Number of 4-week periods (approximately months) from transplant",
    "to measurement. Centred on 1."
  ),

  "weeks4_since_tx_baseline_centered", exp_beta("{XC}"), "$X_C$",
  paste(
    "Expected fold-titre increase for either group at visits 2, 3 and 4",
    "for a 4-week increase in time from transplant.",
    adjusted("weeks4_since_tx_baseline_centered")
  ),
  paste(
    "Number of 4-week periods (approximately months) from transplant",
    "to first visit. Centred on 1."
  ),

  "logtitre_baseline_centered", exp_beta("{BC}"), "$B_C$",
  paste(
    "Expected fold-titre increase for either group at visits 2, 3 and 4",
    "for a 2-fold increase in the baseline titre.",
    adjusted("logtitre_baseline_centered")
  ),
  paste(
    "Log2-baseline titre (i.e. titre at visit 1) shifted to midpoint.",
    "Centred on $\\text{log}_2(5)$."
  ),

  "sd_(Intercept).id", "$r_{\\text{Random}}$", "",
  "A measure of between-subject variability.", "",

  "sd_Observation.Residual", "$r_{\\text{Residual}}$", "",
  "A measure of within-subject variability.", "",
)

write_csv(
  left_join(fits, fits_ref, by = "term"),
  file.path(fit_dir, "fits.csv")
)

fits_ref_rel <- fits_ref %>%
  filter(term %in% fits$term)

fits_ref_vars <- fits_ref_rel %>%
  filter(var_lbl != "") %>%
  mutate(
    varline = paste0(
      var_lbl, " --- ", var_int, ifelse(
        row_number() == max(row_number()), "", "\\\\"
      )
    )
  )

# Variable interpretation
writeLines(
  fits_ref_vars$varline,
  file.path(fit_dir, "vars.tex")
)

# Formula
with(fits_ref_vars, {
  terms <- clean_label(term_lbl)
  vars <- clean_label(var_lbl)
  paste(terms[-1], vars[-1]) %>%
    paste(collapse = " + ") %>%
    paste(vars[1], " = \\beta_0 + ", ., "+ s + e") %>%
    write(file.path(fit_dir, "formula.tex"))
})

# Parameter interpretation
fits_ref_rel %>%
  select(Term = term_lbl, Interpretation = term_int) %>%
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
  column_spec(2, width = "35em") %>%
  write(file.path(fit_dir, "fit-interpret.tex"))
