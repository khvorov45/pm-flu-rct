cat("Extract data for analysis\n")

suppressPackageStartupMessages(library(tidyverse))

data_raw_dir <- here::here("data-raw")
data_dir <- here::here("data")

# Functions ===================================================================

virus_titre_names <- function(...) {
  map(list(...), ~ paste(.x, c("t1", "t2", "t3", "t4"), sep = "_")) %>%
    unlist()
}

save_csv <- function(data, name) {
  write_csv(data, file.path(data_dir, glue::glue("{name}.csv")))
}

# Script ======================================================================

# HI titres

hi <- readxl::read_excel(
  file.path(data_raw_dir, "HI results_HSCT_sent.xlsx"),
  sheet = "Sheet1",
  range = "A12:S79",
  col_names = c(
    "id", "temp", "temp2",
    virus_titre_names("A(H1N1)pdm09", "A(H3N2)", "B Yam", "B Vic")
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
    titre = if_else(titre == "<10", "5", titre) %>% as.integer(),
    logtitre = log(titre),
    logtitre_mid = if_else(titre == 5L, log(5), log(titre) + log(2) / 2),
  ) %>%
  select(-virus_timepoint)

# Subject dates

redcap <- readxl::read_excel(
  file.path(data_raw_dir, "flu raw data export from redcap 9-4-20.xlsx"),
  sheet = 1
) %>%
  mutate(id = str_pad(studyid, 3, pad = "0"))

subjects <- redcap %>%
  filter(
    redcap_event_name == "baseline_v1_arm_1",
    is.na(redcap_repeat_instrument)
  ) %>%
  mutate(
    dob = lubridate::as_date(dob),
    date_x = lubridate::as_date(time_tx),
    date_t1 = lubridate::as_date(vacc_date_1),
    date_t2 = lubridate::as_date(vacc_date_2),
    date_t3 = lubridate::as_date(date_visit_3),
    date_t4 = lubridate::as_date(date_visit_4),
    myeloma = as.integer(cancer_type == 1),
    vac_in_prior_year = as.integer(
      lubridate::year(lubridate::as_date(prior_vacc_date)) >= 2018
    ) %>% replace_na(0L),
  ) %>%
  select(
    id, dob, date_x, contains("date_t"), myeloma, vac_in_prior_year,
    current_therapy
  ) %>%
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
  summarise(mean_offset = mean(days_offset, na.rm = TRUE), .groups = "drop")

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

# Infection events

all_inf_ids <- redcap %>%
  filter(!is.na(redcap_repeat_instrument)) %>%
  pull(id) %>%
  unique()

# Missing infection status over the whole follow-up period because not followed
# for the whole period
mis_inf_ids <- hi %>%
  group_by(id) %>%
  summarise(mis = any(is.na(titre)), .groups = "drop") %>%
  filter(mis) %>%
  pull(id)

subjects_final <- subjects_imp_date %>%
  mutate(
    ili = as.integer(id %in% all_inf_ids),
    ili = if_else(id %in% mis_inf_ids, NA_integer_, ili),
    age_years = (date - dob) / lubridate::dyears(1),
    age_years_centered = age_years - 50,
    weeks4_since_tx = (date - date_x) / lubridate::dweeks(1),
    weeks4_since_tx_centered = weeks4_since_tx - 1,
    timepoint_lbl = factor(
      timepoint,
      levels = 1:4,
      labels = c(
        "Visit 1 (pre-vac1)", "Visit 2 (pre-vac2)", "Visit 3", "Visit 4"
      )
    ),
  ) %>%
  group_by(id) %>%
  mutate(
    days_since_t1 = (date - date[timepoint == 1L]) / lubridate::ddays(1)
  ) %>%
  ungroup() %>%
  select(-dob, -date_x, -date)

# Group assignment
groups <- readxl::read_excel(
  file.path(data_raw_dir, "list of patients.xlsx"),
  sheet = 1,
  range = cellranger::cell_cols(c(1, 2))
) %>%
  mutate(group = factor(
    Arm,
    levels = c("STD", "HD"),
    labels = c("Standard Dose", "High Dose")
  )) %>%
  select(id = SN, group)

all_data <- hi %>%
  inner_join(subjects_final, c("id", "timepoint")) %>%
  inner_join(groups, "id")

# See if there are any duplicates
dups <- all_data %>%
  group_by(id, virus) %>%
  filter(!(all(timepoint == c(1L, 2L, 3L, 4L))))

# Should be 001-068
ids <- all_data$id %>% unique()

# Attach some variables

all_data_extra <- all_data %>%
  group_by(id, virus) %>%
  mutate(
    logtitre_baseline = log2(exp(logtitre[timepoint == 1L])),
    logtitre_baseline_centered = logtitre_baseline - log2(5),
    age_years_baseline_centered = age_years_centered[timepoint == 1L],
    weeks4_since_tx_baseline_centered =
      weeks4_since_tx_centered[timepoint == 1L],
  ) %>%
  ungroup()

save_csv(all_data_extra, "data")

# Modified data

data_titre <- all_data_extra %>%
  filter(timepoint != 1L)
save_csv(data_titre, "titre")

data_ili <- all_data_extra %>%
  filter(timepoint == 1L) %>%
  pivot_wider(names_from = "virus", values_from = contains("titre"))
save_csv(data_ili, "ili")

data_prot_and_conv <- all_data_extra %>%
  group_by(
    virus, id, group, myeloma, vac_in_prior_year, current_therapy,
    age_years_baseline_centered, weeks4_since_tx_baseline_centered
  ) %>%
  summarise(
    seroprotection_before = titre[timepoint == 1L] >= 40L,
    seroprotection = as.integer(titre[timepoint == 3L] >= 40L),
    seroconversion = as.integer(
      titre[timepoint == 3L] / titre[timepoint == 1L] >= 4
    ),
    ratio_v2_v1 = titre[timepoint == 2L] / titre[timepoint == 1L],
    ratio_v3_v1 = titre[timepoint == 3L] / titre[timepoint == 1L],
    ratio_v4_v1 = titre[timepoint == 4L] / titre[timepoint == 1L],
    ratio_v3_v2 = titre[timepoint == 3L] / titre[timepoint == 2L],
    .groups = "drop"
  )
data_seroprotection <- data_prot_and_conv %>%
  filter(!seroprotection_before)
save_csv(data_seroprotection, "seroprotection")
save_csv(data_prot_and_conv, "seroconversion")

data_conv_and_prot_combined <- all_data_extra %>%
  filter(virus != "B Vic") %>%
  group_by(
    id, group, myeloma, vac_in_prior_year, current_therapy,
    age_years_baseline_centered, weeks4_since_tx_baseline_centered
  ) %>%
  summarise(
    n_seroprotection_before = sum(titre[timepoint == 1L] >= 40L),
    n_seroprotection = sum(titre[timepoint == 3L] >= 40L),
    n_seroconversion = sum(
      titre[timepoint == 3L] / titre[timepoint == 1L] >= 4
    ),
    .groups = "drop"
  ) %>%
  map_dfr(1:3, function(n_prot, df) {
    df %>%
      mutate(
        seroprotection_before = n_seroprotection_before >= n_prot,
        seroprotection = as.integer(n_seroprotection >= n_prot),
        seroconversion = as.integer(n_seroconversion >= n_prot),
        n_prot = n_prot
      )
  }, .)
data_seroprotection_combined <- data_conv_and_prot_combined %>%
  filter(!seroprotection_before)
save_csv(data_seroprotection_combined, "seroprotection_combined")
save_csv(data_conv_and_prot_combined, "seroconversion_combined")

# Make sure we've got 1 row per individual
stopifnot(all(data_ili$id == unique(data_titre$id)))
data_seroprotection_test <- data_seroprotection %>%
  group_by(virus) %>%
  summarise(
    ids = length(id),
    unique_ids = length(unique(id)),
    .groups = "drop"
  ) %>%
  filter(ids != unique_ids)
stopifnot(nrow(data_seroprotection_test) == 0)
data_seroprotection_combined_test <- data_seroprotection_combined %>%
  group_by(n_prot) %>%
  summarise(
    ids = length(id),
    unique_ids = length(unique(id)),
    .groups = "drop"
  ) %>%
  filter(ids != unique_ids)
stopifnot(nrow(data_seroprotection_combined_test) == 0)

# Adverse events

# The events
adverse_events <- redcap %>%
  select(id, matches("vacc_reaction_[1|2]")) %>%
  pivot_longer(
    -id,
    names_to = "adverse_event", values_to = "adverse_outcome"
  ) %>%
  mutate(
    vaccine_index = str_replace(
      adverse_event, "vacc_reaction_(\\d)___\\d", "\\1"
    ) %>% as.integer(),
    adverse_event = str_replace(adverse_event, ".*(\\d)$", "\\1") %>%
      recode(
        "1" = "None", "2" = "Myalgia", "3" = "Swelling", "4" = "Erythema",
        "5" = "Malaise", "6" = "Headache", "7" = "Fever"
      )
  ) %>%
  filter(adverse_outcome == 1, adverse_event != "None") %>%
  select(id, vaccine_index, adverse_event)

# Their severity
adverse_severity <- redcap %>%
  select(id, contains("severity")) %>%
  select(-contains("disease"), -diabetes_severity, -severity_ctcae) %>%
  pivot_longer(-id, names_to = "adverse_event", values_to = "severity") %>%
  filter(!is.na(severity)) %>%
  mutate(
    vaccine_index = str_replace(adverse_event, ".*(\\d)$", "\\1") %>%
      as.integer(),
    adverse_event = str_replace(adverse_event, "_severity_\\d$", "") %>%
      tools::toTitleCase(),
    severity = recode(severity, "1" = "Mild", "2" = "Moderate", "3" = "Severe")
  )

# Make sure events match
stopifnot(identical(
  setdiff(adverse_events$adverse_event, adverse_severity$adverse_event),
  character(0)
))
stopifnot(identical(
  setdiff(adverse_severity$adverse_event, adverse_events$adverse_event),
  character(0)
))

adverse_events_full <- full_join(
  adverse_events, adverse_severity,
  by = c("id", "vaccine_index", "adverse_event")
)

# Make sure all events have associated severity
stopifnot(nrow(filter_all(adverse_events_full, is.na)) == 0L)

adverse_events_full_groups <- inner_join(
  adverse_events_full, groups,
  by = "id"
)
stopifnot(nrow(adverse_events_full_groups) == nrow(adverse_events_full))

save_csv(adverse_events_full_groups, "adverse_events")
