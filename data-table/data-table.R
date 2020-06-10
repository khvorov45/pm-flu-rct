cat("Make data table")

library(tidyverse)
library(kableExtra)

data_dir <- here::here("data")
data_table_dir <- here::here("data-table")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

save_csv <- function(table, name) {
  write_csv(table, file.path(data_table_dir, paste0(name, ".csv")))
  table
}

save_table <- function(table, name) {
  write(table, file.path(data_table_dir, paste0(name, ".tex")))
}

# Script ======================================================================

data <- read_data()

miss_counts <- data %>%
  group_by(virus, timepoint_lbl) %>%
  summarise(
    n_nomiss = sum(!is.na(titre))
  )

miss_counts_tbl <- miss_counts %>%
  pivot_wider(names_from = "virus", values_from = n_nomiss) %>%
  rename(Timepoint = timepoint_lbl) %>%
  save_csv("nobs") %>%
  kable(
    format = "latex",
    caption =
      "Number of observations at the four time points for the four viruses.",
    label = "nobs",
    booktabs = TRUE,
    align = "lcccc"
  ) %>%
  kable_styling(latex_options = "striped") %>%
  save_table("nobs")
