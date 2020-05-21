read_data <- function() {
  read_csv(
    file.path("data/data.csv"),
    col_types = cols(
      titre = col_integer(),
      timepoint = col_integer(),
      date_imputed = col_integer(),
      timepoint_lbl = col_factor(
        c("Pre-V1", "Pre-V2", "Post-V2 Visit 1", "Post-V2 Visit 2")
      )
    )
  )
}
