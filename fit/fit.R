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
      #+ myeloma
      #+ current_therapy
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
    ili ~ group
    # myeloma
    #+ current_therapy
    + age_years_baseline_centered,
    #+ weeks4_since_tx_baseline_centered,
    binomial,
    data
  ) %>%
    broom::tidy()
}

fit_model_seroprotection <- function(data) {
  glm(
    seroprotection ~ group
    #+ myeloma
    #+ current_therapy
    + age_years_baseline_centered,
    #+ weeks4_since_tx_baseline_centered,
    binomial,
    data
  ) %>%
    broom::tidy()
}

fit_model_seroprotection_combined <- function(data, key) {
  formulas <- list(
    seroprotection ~ group
    #+ myeloma
    #+ current_therapy
    + age_years_baseline_centered,
    #+ weeks4_since_tx_baseline_centered,
    seroprotection ~ group
    #+ myeloma
    #+ current_therapy
    + age_years_baseline_centered,
    #+ weeks4_since_tx_baseline_centered,
    seroprotection ~ group
    #+ myeloma
    #+ current_therapy
    + age_years_baseline_centered
    #+ weeks4_since_tx_baseline_centered,
  )
  glm(
    formulas[[key$n_prot]],
    binomial,
    data
  ) %>%
    broom::tidy()
}

fit_model_seroconversion <- function(data) {
  glm(
    seroconversion ~ group
    #+ myeloma
    #+ current_therapy
    + age_years_baseline_centered,
    #+ weeks4_since_tx_baseline_centered,
    binomial,
    data
  ) %>%
    broom::tidy()
}

fit_model_seroconversion_combined <- function(data, key) {
  glm(
    seroconversion ~ group
    #+ myeloma
    + vac_in_2018
      + vac_in_2019
      + current_therapy
      + age_years_baseline_centered
      + weeks4_since_tx_baseline_centered,
    binomial,
    data
  ) %>%
    broom::tidy()
}

fit_seroconversion_individual <- function(data, vars) {
  vars <- c(
    "group",
    "vac_in_2018",
    "vac_in_2019",
    "current_therapy",
    "age_years_baseline_centered",
    "weeks4_since_tx_baseline_centered"
  )
  map_dfr(
    vars,
    ~ broom::tidy(glm(
      as.formula(glue::glue("seroconversion ~ {.x}")),
      binomial,
      data
    )) %>%
      mutate(variable_name = .x)
  )
}

summ_seroconversion_individual <- function(data) {
  p <- function(n, c) {
    prop <- signif(n / c * 100, 2)
    glue::glue("{n} ({prop}\\%)")
  }
  m <- function(vec) {
    glue::glue("{signif(mean(vec), 2)} $\\pm$ {signif(sd(vec), 2)}")
  }
  data %>%
    filter(!is.na(seroconversion)) %>%
    group_by(seroconversion) %>%
    summarise(
      group = p(sum(group == "High Dose"), n()),
      vac_in_2018 = p(sum(vac_in_2018), n()),
      vac_in_2019 = p(sum(vac_in_2019), n()),
      current_therapy = p(sum(current_therapy), n()),
      age_years_baseline_centered = m(age_years_baseline_centered + 50),
      weeks4_since_tx_baseline_centered = m(weeks4_since_tx_baseline_centered),
      ntot = n(),
      .groups = "keep"
    ) %>%
    mutate(
      seroconversion = recode(
        seroconversion,
        "1" = glue::glue("Seroconverted (n={ntot})"),
        "0" = glue::glue("Did not seroconvert (n={ntot})")
      )
    ) %>%
    ungroup() %>%
    select(-ntot) %>%
    pivot_longer(
      c(-seroconversion),
      names_to = "variable_name", values_to = "data_summ"
    ) %>%
    pivot_wider(names_from = "seroconversion", values_from = "data_summ")
}

gen_b0_int <- function(fits, what, when) {
  all <- list(
    "myeloma" = "cancer other than myeloma",
    "vac_in_2018" = "not vaccinated in 2018",
    "vac_in_2019" = "not vaccinated in 2019",
    "current_therapy" = "not receiving therapy",
    "age_years_centered" = "age 50",
    "age_years_baseline_centered" = "age 50",
    "weeks4_since_tx_centered" = "1 week since transplant",
    "weeks4_since_tx_baseline_centered" = "1 week since transplant",
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
    "vac_in_2018" = "vaccination in 2018",
    "vac_in_2019" = "vaccination in 2019",
    "current_therapy" = "current therapy status",
    "age_years_centered" = "age",
    "age_years_baseline_centered" = "age",
    "weeks4_since_tx_centered" = "time from transplant",
    "weeks4_since_tx_baseline_centered" = "time from transplant",
    "logtitre_baseline_centered" = "baseline titre"
  )
  adj <- adj[names(adj) %in% fits$term]
  if (!is.null(this)) adj <- adj[!names(adj) %in% this]
  if (length(adj) == 0) {
    return("")
  }
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

    "vac_in_2018", fun_beta("{PV2018}", beta_fun_name), "$P_{2018}$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for subjects vaccinated in 2018 ",
      "compared to subjects not vaccinated in 2018. ",
      adjusted(fits, "vac_in_2018")
    ),
    paste(
      "Indicator of whether the subject was vaccinated in 2018",
      "(1) or not (0)."
    ),

    "vac_in_2019", fun_beta("{PV2019}", beta_fun_name), "$P_{2019}$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for subjects vaccinated in 2019 ",
      "compared to subjects not vaccinated in 2019. ",
      adjusted(fits, "vac_in_2019")
    ),
    paste(
      "Indicator of whether the subject was vaccinated in 2019",
      "(1) or not (0)."
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
      "for a 1-week increase in time from transplant. ",
      adjusted(fits, "weeks4_since_tx_centered")
    ),
    paste(
      "Number of weeks from transplant ",
      "to measurement. Centred on 1."
    ),

    "weeks4_since_tx_baseline_centered",
    fun_beta("{XC}", beta_fun_name), "$X_C$",
    glue::glue(
      "Expected {what_other} for either group{when_other} ",
      "for a 1-week increase in time from transplant. ",
      adjusted(fits, "weeks4_since_tx_baseline_centered")
    ),
    paste(
      "Number of weeks from transplant",
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

save_csv <- function(data, name) {
  write_csv(data, file.path(fit_dir, glue::glue("{name}.csv")))
}

save_fits <- function(fits, fits_ref, name) {
  write_csv(
    left_join(fits, fits_ref, by = "term"),
    file.path(fit_dir, glue::glue("fits-{name}.csv"))
  )
}

save_table <- function(table, name) {
  write(table, file.path(fit_dir, paste0(name, ".tex")))
}

# Script ======================================================================

# Prepare data

data_titre <- read_data("titre")
data_ili <- read_data("ili")
data_seroprotection <- read_data("seroprotection")
data_seroprotection_combined <- read_data("seroprotection_combined")
data_seroconversion <- read_data("seroconversion")
data_seroconversion_combined <- read_data("seroconversion_combined")

# Fit models

fits_titre <- data_titre %>%
  group_by(virus) %>%
  group_modify(~ fit_model(.x)) %>%
  mutate(p.value = 2 * pnorm(-abs(statistic)))

fits_ili <- data_ili %>%
  fit_model_ili()

fits_seroprotection <- data_seroprotection %>%
  group_by(virus) %>%
  group_modify(~ fit_model_seroprotection(.x))

fits_seroprotection_combined <- data_seroprotection_combined %>%
  group_by(n_prot) %>%
  group_modify(fit_model_seroprotection_combined)

fits_seroconversion <- data_seroconversion %>%
  group_by(virus) %>%
  group_modify(~ fit_model_seroconversion(.x))

fits_seroconversion_combined <- data_seroconversion_combined %>%
  group_by(n_prot) %>%
  group_modify(fit_model_seroconversion_combined)

fit_seroconversion_combined3_individual <- data_seroconversion_combined %>%
  filter(n_prot == 3) %>%
  fit_seroconversion_individual() %>%
  filter(term != "(Intercept)") %>%
  inner_join(
    summ_seroconversion_individual(
      filter(data_seroconversion_combined, n_prot == 3)
    ),
    by = "variable_name"
  )

# Reference tables

fits_ref_titre <- gen_fits_ref(
  fits_titre, "titre", "fold-titre change", " at visit 2", " at any visit",
  "\\text{exp}", "$T_{l,m}$",
  "Natural log-titre shifted to the midpoint."
)
fits_ref_ili <- gen_fits_ref(
  fits_ili, "odds of infection", "odds ratio of infection", "", "",
  "\\text{exp}", "$I$",
  "Indicator of infection (1) or no infection (0) during the study period."
)
fits_ref_seroprotection <- gen_fits_ref(
  fits_seroprotection, "odds of achieving seroprotection",
  "odds ratio of achieving seroprotection", "", "",
  "\\text{exp}", "$S$",
  "Indicator of seroprotection (1) or no seroprotection (0) at visit 3."
)
fits_ref_seroprotection_combined <- gen_fits_ref(
  fits_seroprotection_combined,
  paste(
    "odds of achieving seroprotection",
    "against at least the given number of antigens"
  ),
  paste(
    "odds ratio of achieving seroprotection",
    "against at least the given number of antigens"
  ), "", "",
  "\\text{exp}", "$S_C$",
  paste(
    "Indicator of seroprotection against at least",
    "the given number of antigens (1) or not (0) at visit 3."
  )
)
fits_ref_seroconversion <- gen_fits_ref(
  fits_seroconversion, "odds of achieving seroconversion",
  "odds ratio of achieving seroconversion", "", "",
  "\\text{exp}", "$S$",
  "Indicator of seroconversion (1) or no seroconversion (0) at visit 3."
)
fits_ref_seroconversion_combined <- gen_fits_ref(
  fits_seroconversion_combined,
  paste(
    "odds of achieving seroconversion",
    "against at least the given number of antigens"
  ),
  paste(
    "odds ratio of achieving seroconversion",
    "against at least the given number of antigens"
  ), "", "",
  "\\text{exp}", "$S_C$",
  paste(
    "Indicator of seroconversion against at least",
    "the given number of antigens (1) or not (0) at visit 3."
  )
)

save_fits(fits_titre, fits_ref_titre, "titre")
save_fits(fits_ili, fits_ref_ili, "ili")
save_fits(fits_seroprotection, fits_ref_seroprotection, "seroprotection")
save_fits(
  fits_seroprotection_combined,
  fits_ref_seroprotection_combined, "seroprotection_combined"
)
save_fits(fits_seroconversion, fits_ref_seroconversion, "seroconversion")
save_fits(
  fits_seroconversion_combined,
  fits_ref_seroconversion_combined, "seroconversion_combined"
)
save_fits(
  fit_seroconversion_combined3_individual,
  fits_ref_seroconversion_combined, "seroconversion_combined3_individual"
)

fits_ref_all <- list(
  "titre" = fits_ref_titre %>% filter(term %in% fits_titre$term),
  "ili" = fits_ref_ili %>% filter(term %in% fits_ili$term),
  "seroprotection" = fits_ref_seroprotection %>%
    filter(term %in% fits_seroprotection$term),
  "seroprotection_combined" = fits_ref_seroprotection_combined %>%
    filter(term %in% fits_seroprotection_combined$term),
  "seroconversion" = fits_ref_seroconversion %>%
    filter(term %in% fits_seroconversion$term),
  "seroconversion_combined" = fits_ref_seroconversion_combined %>%
    filter(term %in% fits_seroconversion_combined$term)
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
    save_csv(fits_ref_vars, glue::glue("{name}-interpretation"))
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
iwalk(
  fits_ref_all_vars, function(fits_ref_vars, name) {
    if (name == "titre") {
      return()
    }
    with(fits_ref_vars, {
      terms <- clean_label(term_lbl)
      vars <- clean_label(var_lbl)
      paste(terms[-1], vars[-1]) %>%
        paste(collapse = " + ") %>%
        paste(
          paste0("E(\\text{logit}(P(", vars[1], "=1)))"), " = \\beta_0 + ", .
        ) %>%
        write(file.path(fit_dir, glue::glue("formula-{name}.tex")))
    })
  }
)

# Parameter interpretation
iwalk(fits_ref_all, function(fits_ref_rel, name) {
  models <- c(
    "titre" = "titre",
    "ili" = "ILI",
    "seroprotection" = "seroprotection",
    "seroprotection_combined" = "combined seroprotection",
    "seroconversion" = "seroconversion",
    "seroconversion_combined" = "combined seroconversion"
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
