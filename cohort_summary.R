# script to explore the cohort qualities
pacman::p_load(tidyverse, mgcv, gratia, glue, modelsummary)

source(here::here("utils.R"))
source(here::here("load_data.R"))

ggplot2::theme_set(ggplot2::theme_bw())


output_location <- fs::dir_create(here::here(
  "outputs", "tables", "cohort"
))

study_start_date <- "2023-11-14"
study_end_date <- "2024-03-07"


outcomes <- c(
  "symptoms_new_cough",
  "symptoms_new_fever",
  "symptoms_new_sore_throat",
  "symptoms_new_short_breath",
  "symptoms_new_loss_smell"
)

# we want to show critical information about the cohort, such as demographics,
# symptom onset information and how many waves they completed.
cohort <- symptoms |>
  dplyr::select(
    submitted_at, participant_id, sex, age_group, rgn21nm, imd_dec, dplyr::all_of(outcomes)
  ) |>
  dplyr::mutate(imd_quint = factor(ceiling(imd_dec / 2)), .keep = "unused")

cohort |>
  dplyr::select(-c(dplyr::any_of(outcomes), participant_id)) |>
  dplyr::mutate(age_group = stringr::str_remove(age_group, " years")) |>
  dplyr::rename(
    "Sex" = sex,
    "Age group" = age_group,
    "Region" = rgn21nm,
    "IMD quintile" = imd_quint
  ) |>
  modelsummary::datasummary_skim(type = "categorical", output = fs::path(output_location, "demographics.png"))

cohort |>
  dplyr::select(dplyr::any_of(outcomes)) |>
  dplyr::rename_with(.fn = \(x) {
    stringr::str_remove(x, "symptoms_new_") |>
      stringr::str_replace("_", " ") |>
      stringr::str_to_sentence()
  }) |>
  modelsummary::datasummary_skim(type = "categorical", output = fs::path(output_location, "symptoms.png"))

# there are 4 possible waves over the study period
# we want to know how many each participant responded to and summarise
# which indicates how independently we can assume the samples are
cohort |>
  dplyr::summarise(
    waves = factor(dplyr::n(), levels = c(1, 2, 3, 4)),
    .by = "participant_id"
  ) |>
  modelsummary::datasummary_skim(type = "categorical", output = fs::path(output_location, "waves.png"))

# how many unique participants
dplyr::n_distinct(cohort$participant_id)

# how many submissions
nrow(cohort)
