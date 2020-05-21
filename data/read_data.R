read_data <- function() {
  read_csv(
    file.path("data/data.csv"),
    col_types = cols(
      titre = col_integer(),
      timepoint = col_integer(),
      date_imputed = col_integer()
    )
  )
}
