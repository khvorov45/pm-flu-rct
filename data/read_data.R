read_data <- function() {
  read_csv(
    file.path("data/data.csv"),
    col_types = cols(
      titre = col_integer(),
      timepoint = col_integer(),
      date_imputed = col_integer(),
      timepoint_lbl = col_factor(
        c("Visit 1 (pre-vac1)", "Visit 2 (pre-vac2)", "Visit 3", "Visit 4")
      ),
      group = col_factor(c("Standard Dose", "High Dose"))
    )
  )
}
