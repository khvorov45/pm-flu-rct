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

# HI titres

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

# Subject dates

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
  mutate(timepoint = str_replace(timepoint, "date_t", "") %>% as.integer())

# Impute dates

# We are missing T4 date mostly
subjects_miss_date <- subjects %>% filter(is.na(date))

# But some subjects actually had blood taken at timepoint 4
t_miss_but_blood_taken <- inner_join(
  subjects_miss_date, hi, c("id", "timepoint")
) %>% filter(!is.na(titre))

mean_times <- subjects %>%
  group_by(id) %>%
  mutate(days_offset = date - date[timepoint == 1L]) %>%
  group_by(timepoint) %>%
  summarise(mean_offset = mean(days_offset, na.rm = TRUE))

subjects_imp_date <- subjects %>%
  inner_join(mean_times, "timepoint") %>%
  group_by(id) %>%
  mutate(
    date_imputed = is.na(date) %>% as.integer(),
    date = if_else(
      is.na(date), date[timepoint == 1L] + mean_offset, date
    )
  ) %>%
  ungroup() %>%
  select(-mean_offset)

# Check for wrong imputation
subjects_imp_date %>%
  group_by(id) %>%
  filter(any(date < date[timepoint == 1L]))

subjects_final <- subjects_imp_date %>%
  mutate(
    age_years = (date - dob) / lubridate::dyears(1),
    days_since_tx = (date - date_x) / lubridate::ddays(1),
    timepoint_lbl = factor(
      timepoint,
      levels = 1:4,
      labels = c("Pre-V1", "Pre-V2", "Post-V2 Visit 1", "Post-V2 Visit 2")
    ),
  ) %>%
  group_by(id) %>%
  mutate(
    days_since_t1 = (date - date[timepoint == 1L]) / lubridate::ddays(1)
  ) %>%
  ungroup() %>%
  select(-dob, -date_x, -date)

groups <- readxl::read_excel(
  file.path(data_raw_dir, "list of patients.xlsx"),
  sheet = 1,
  range = cellranger::cell_cols(c(1, 2))
) %>%
  mutate(group = recode(Arm, "HD" = "High Dose", "STD" = "Standard Dose")) %>%
  select(id = SN, group)

all_data <- hi %>%
  inner_join(subjects_final, c("id", "timepoint")) %>%
  inner_join(groups, "id")

# See if there are any duplicates
all_data %>%
  group_by(id, virus) %>%
  filter(!(all(timepoint == c(1L, 2L, 3L, 4L))))

# Should be 001-068
all_data$id %>% unique()

write_csv(all_data, file.path(data_dir, "data.csv"))
