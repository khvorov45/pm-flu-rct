# Extract the data for analysis

library(tidyverse)

data_raw_dir <- here::here("data-raw")
data_dir <- here::here("data")

# Functions ===================================================================

virus_titre_names <- function(...) {
  map(list(...), ~ paste(.x, c("t1", "t2", "t3", "t4"), sep = "_")) %>%
    unlist()
}

# Script ======================================================================

hi <- readxl::read_excel(
  file.path(data_raw_dir, "HI results_HSCT_sent.xlsx"),
  sheet = "Sheet1",
  range = "A12:S79",
  col_names = c(
    "id", "temp", "temp2",
    virus_titre_names("A Michigan", "A Brisbane", "B Yam", "B Vic")
  ),
  col_types = "text",
  na = c("", "No sample")
) %>%
  select(-contains("temp")) %>%
  pivot_longer(-id, names_to = "virus_timepoint", values_to = "titre") %>%
  mutate(
    timepoint = str_replace(virus_timepoint, "^.*(\\d)$", "\\1") %>%
      as.integer(),
    virus = str_replace(virus_timepoint, "(^.*)_t\\d$", "\\1"),
    titre = if_else(titre == "<10", "5", titre) %>% as.integer()
  ) %>%
  select(-virus_timepoint)

subjects <- readxl::read_excel(
  file.path(data_raw_dir, "flu raw data export from redcap 9-4-20.xlsx"),
  sheet = 1
) %>%
  filter(
    redcap_event_name == "baseline_v1_arm_1",
    is.na(redcap_repeat_instrument)
  ) %>%
  mutate(
    id = str_pad(studyid, 3, pad = "0"),
    dob = lubridate::as_date(dob),
    date_x = lubridate::as_date(time_tx),
    date_t1 = lubridate::as_date(vacc_date_1),
    date_t2 = lubridate::as_date(vacc_date_2),
    date_t3 = lubridate::as_date(date_visit_3),
    date_t4 = lubridate::as_date(date_visit_4),
  ) %>%
  select(id, dob, date_x, contains("date_t")) %>%
  pivot_longer(
    contains("date_t"),
    names_to = "timepoint", values_to = "date"
  ) %>%
  mutate(
    timepoint = str_replace(timepoint, "date_t", "") %>% as.integer(),
    age_years = (date - dob) / lubridate::dyears(1),
    days_since_tx = (date - date_x) / lubridate::ddays(1),
  ) %>%
  select(-dob, -date_x)

groups <- readxl::read_excel(
  file.path(data_raw_dir, "list of patients.xlsx"),
  sheet = 1,
  range = cellranger::cell_cols(c(1, 2))
) %>%
  mutate(group = recode(Arm, "HD" = "High Dose", "STD" = "Standard Dose")) %>%
  select(id = SN, group)

all_data <- hi %>%
  inner_join(subjects, c("id", "timepoint")) %>%
  inner_join(groups, "id")

# See if there are any duplicates
all_data %>%
  group_by(id, virus) %>%
  filter(n() != 4)

# Should be 001-068
all_data$id %>% unique()

write_csv(all_data, file.path(data_dir, "data.csv"))
