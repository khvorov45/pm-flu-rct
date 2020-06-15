cat("Fit a model\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(kableExtra)
})

data_dir <- here::here("data")
fit_dir <- here::here("fit")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

fun_beta <- function(beta_name, fun_name) {
  if (str_detect(beta_name, "\\{")) {
    beta_name <- str_replace(beta_name, "\\{(.*)\\}", "\\{\\\\text\\{\\1\\}\\}")
  }
  paste0("$", fun_name, "(\\beta_", beta_name, ")$")
}

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

fit_model_ili <- function(data) {
  glm(
    ili ~ group +
      myeloma
      + vac_in_prior_year
      + current_therapy
      + age_years_baseline_centered
      + weeks4_since_tx_baseline_centered,
    binomial,
    data
  ) %>%
    broom::tidy()
}

gen_b0_int <- function(fits, what, when) {
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
  glue::glue(
    "Expected {what} for the standard dose group{when}, ",
    paste(all, collapse = ", "), "."
  )
}

adjusted <- function(fits, this = NULL) {
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

gen_fits_ref <- function(fits, what_int, what_other, when_int, when_other,
                         beta_fun_name, var_lbl_int, var_int_int) {
  all <- tribble(
    ~term, ~term_lbl, ~var_lbl, ~term_int, ~var_int,

    "(Intercept)", fun_beta("0", beta_fun_name), var_lbl_int,
    gen_b0_int(fits, what_int, when_int), var_int_int,

    "groupHigh Dose", fun_beta("{HD}", beta_fun_name), "$G_H$",
    glue::glue(
      "Expected {what_other} for the high dose group",
      " as compared to the standard dose group",
      "{when_other}. ",
      adjusted(fits)
    ),
    paste(
      "Indicator of whether the observation belongs to the",
      "high-dose group (1) or the standard-dose group (0)."
    ),

    "timepoint_lblVisit 3", fun_beta("{V3}", beta_fun_name), "$V_3$",
    glue::glue(
      "Expected {what_other} for either group",
      " at visit 3 as compared to visit 2. ",
      adjusted(fits)
    ),
    paste(
      "Indicator of whether the observation belongs to",
      "visit 3 (1) or not (0)."
    ),

    "timepoint_lblVisit 4", fun_beta("{V4}", beta_fun_name), "$V_4$",
    glue::glue(
      "Expected {what_other} for either group",
      " at visit 4 as compared to visit 2. ",
      adjusted(fits)
    ),
    paste(
      "Indicator of whether the observation belongs to",
      "visit 4 (1) or not (0)."
    ),

    "myeloma", fun_beta("M", beta_fun_name), "$M$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for subjects with myeloma compared to subjects with other cancer type. ",
      adjusted(fits, "myeloma")
    ),
    "Indicator of whether the subject has myeloma (1) or not (0)",

    "vac_in_prior_year", fun_beta("{PV}", beta_fun_name), "$P$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for subjects last vaccinated no earlier than 2018 ",
      "compared to subjects last vaccinated earlier ",
      "than 2018 or never vaccinated. ",
      adjusted(fits, "vac_in_prior_year")
    ),
    paste(
      "Indicator of whether the subject was vaccinated no earlier than 2018",
      "(1) or either vaccinated earlier than 2018 or never vaccinated (0)."
    ),

    "current_therapy", fun_beta("{R}", beta_fun_name), "$R$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for subjects receiving therapy as compared to subjects not ",
      "receiving therapy. ",
      adjusted(fits, "current_therapy")
    ),
    paste(
      "Indicator of whether the subject was receiving therapy (1)",
      "or not (0)."
    ),

    "age_years_centered", fun_beta("{AC}", beta_fun_name), "$A_C$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for 1 year increase in age. ",
      adjusted(fits, "age_years_centered")
    ),
    "Age in years at first visit. Centred on 50.",

    "age_years_baseline_centered", fun_beta("{AC}", beta_fun_name), "$A_C$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for 1 year increase in age. ",
      adjusted(fits, "age_years_baseline_centered")
    ),
    "Age in years at time of measurement. Centred on 50.",

    "weeks4_since_tx_centered", fun_beta("{XC}", beta_fun_name), "$X_C$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for a 4-week increase in time from transplant. ",
      adjusted(fits, "weeks4_since_tx_centered")
    ),
    paste(
      "Number of 4-week periods (approximately months) from transplant ",
      "to measurement. Centred on 1."
    ),

    "weeks4_since_tx_baseline_centered",
    fun_beta("{XC}", beta_fun_name), "$X_C$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for a 4-week increase in time from transplant. ",
      adjusted(fits, "weeks4_since_tx_baseline_centered")
    ),
    paste(
      "Number of 4-week periods (approximately months) from transplant",
      "to first visit. Centred on 1."
    ),

    "logtitre_baseline_centered", fun_beta("{BC}", beta_fun_name), "$B_C$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      " for a 2-fold increase in the baseline titre. ",
      adjusted(fits, "logtitre_baseline_centered")
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
  filter(all, term %in% fits$term)
}

clean_label <- function(labels, patt = "\\\\text\\{exp\\}\\((.*)\\)") {
  str_replace_all(labels, "\\$", "") %>%
    str_replace_all(patt, "\\1")
}

# Script ======================================================================

# Prepare data

data <- read_data() %>%
  group_by(id, virus) %>%
  mutate(
    logtitre_baseline = log2(exp(logtitre[timepoint == 1L])),
    logtitre_baseline_centered = logtitre_baseline - log2(5),
    age_years_baseline_centered = age_years_centered[timepoint == 1L],
    weeks4_since_tx_baseline_centered =
      weeks4_since_tx_centered[timepoint == 1L],
  ) %>%
  ungroup()

data_titre <- data %>%
  filter(timepoint != 1L)

data_ili <- data %>%
  filter(timepoint == 1L) %>%
  pivot_wider(names_from = "virus", values_from = contains("titre"))

# Make sure we've got 1 row per individual
stopifnot(all(data_ili$id == unique(data_titre$id)))

# Fit models

fits_titre <- data_titre %>%
  group_by(virus) %>%
  group_modify(~ fit_model(.x))

fits_ili <- data_ili %>%
  fit_model_ili()

# Reference tables

fits_ref_titre <- gen_fits_ref(
  fits_titre, "titre", "fold-titre change", " at visit 2", " at any visit",
  "\\text{exp}", "$T_{l,m}$",
  "Natural log-titre shifted to the midpoint."
)
fits_ref_ili <- gen_fits_ref(
  fits_ili, "odds of infection", "odds ratio of infection", "", "",
  "\\text{exp}", "$I$",
  "Indicator of infection (1) or no infection (0) during the study period"
)

write_csv(
  left_join(fits_titre, fits_ref_titre, by = "term"),
  file.path(fit_dir, "fits-titre.csv")
)

write_csv(
  left_join(fits_ili, fits_ref_ili, by = "term"),
  file.path(fit_dir, "fits-ili.csv")
)

fits_ref_all <- list(
  "titre" = fits_ref_titre %>% filter(term %in% fits_titre$term),
  "ili" = fits_ref_ili %>% filter(term %in% fits_ili$term)
)

fits_ref_all_vars <- map(fits_ref_all, function(fits_ref) {
  fits_ref %>%
    filter(var_lbl != "") %>%
    mutate(
      varline = paste0(
        var_lbl, " --- ", var_int, ifelse(
          row_number() == max(row_number()), "", "\\\\"
        )
      )
    )
})

# Variable interpretation
iwalk(
  fits_ref_all_vars, function(fits_ref_vars, name) {
    writeLines(
      fits_ref_vars$varline,
      file.path(fit_dir, glue::glue("vars-{name}.tex"))
    )
  }
)

# Formula
with(fits_ref_all_vars[["titre"]], {
  terms <- clean_label(term_lbl)
  vars <- clean_label(var_lbl)
  paste(terms[-1], vars[-1]) %>%
    paste(collapse = " + ") %>%
    paste(vars[1], " = \\beta_0 + ", ., "+ s + e") %>%
    write(file.path(fit_dir, "formula-titre.tex"))
})
with(fits_ref_all_vars[["ili"]], {
  terms <- clean_label(term_lbl)
  vars <- clean_label(var_lbl)
  paste(terms[-1], vars[-1]) %>%
    paste(collapse = " + ") %>%
    paste(paste0("E(\\text{logit}P(", vars[1], "=1))"), " = \\beta_0 + ", .) %>%
    write(file.path(fit_dir, "formula-ili.tex"))
})

# Parameter interpretation
iwalk(fits_ref_all, function(fits_ref_rel, name) {
  models <- c(
    "titre" = "titre",
    "ili" = "ILI"
  )
  fits_ref_rel %>%
    select(Term = term_lbl, Interpretation = term_int) %>%
    kable(
      format = "latex",
      caption = glue::glue(
        "Interpretation of the {models[[name]]} model parameters."
      ),
      label = glue::glue("estimates-interpret-{name}"),
      escape = FALSE,
      booktabs = TRUE,
      align = "ll"
    ) %>%
    kable_styling(
      latex_options = "striped"
    ) %>%
    column_spec(2, width = "35em") %>%
    write(file.path(fit_dir, glue::glue("fit-interpret-{name}.tex")))
})
