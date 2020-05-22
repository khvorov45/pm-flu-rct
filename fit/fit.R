# Fit a model

library(tidyverse)

data_dir <- here::here("data")
fit_dir <- here::here("fit")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

# Script ======================================================================

data <- read_data()

lm(logtitre ~ group + timepoint_lbl, data)
