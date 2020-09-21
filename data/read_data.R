read_data <- function(name) {
  suppressWarnings(
    read_csv(
      file.path("data", glue::glue("{name}.csv")),
      col_types = cols(
        titre = col_integer(),
        timepoint = col_integer(),
        date_imputed = col_integer(),
        myeloma = col_integer(),
        vac_in_2018 = col_integer(),
        vac_in_2019 = col_integer(),
        timepoint_lbl = col_factor(
          c("Visit 1 (pre-vac1)", "Visit 2 (pre-vac2)", "Visit 3", "Visit 4")
        ),
        group = col_factor(c("Standard Dose", "High Dose")),
        ili = col_integer(),
        severity = col_factor(c("Mild", "Moderate", "Severe")),
        vaccine_index = col_integer(),
        seroconversion = col_integer()
      )
    )
  )
}
