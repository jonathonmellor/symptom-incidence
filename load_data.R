# script to load and clean all needed WCIS data

pacman::p_load(tidyverse, brms, mgcv, gratia, geomtextpath, sn, glue, ggforce, tidybayes, gratia)

box::use(
  # DATABASE ACCESS - REMOVED
)

print("running load_data script")


source(here::here("utils.R"))

ggplot2::theme_set(ggplot2::theme_bw())

study_start_date <- "2023-11-14"
study_end_date <- "2024-03-07"


## Load data ####

main_survey_path_pattern <- "" # PATH REMOVED FOR CONFIDENTIALITY
main_survey_paths <- file_access$find_latest_file(main_survey_path_pattern)
prepared_data <- s3read_using(
  read_rds,
  object = main_survey_paths
)

main_survey <- prepared_data$main_survey

poststrat <- prepared_data$population

colnames(main_survey)

symptoms <- main_survey |>
  dplyr::select(
    submitted_at, participant_id, sex,
    age, rgn21nm, imd_dec, receive_date,
    dplyr::starts_with("window"),
    dplyr::starts_with("symptoms")
  ) |>
  # relative date variables
  dplyr::mutate(
    window_submission = submitted_at - window_start_date,
    symptoms_submission = submitted_at - symptoms_onset,
    wday = lubridate::wday(submitted_at, label = TRUE)
  ) |>
  # create flag for any symptoms
  dplyr::mutate(
    symptoms_new_any = (symptoms_new_fever | symptoms_new_tiredness | symptoms_new_headache | symptoms_new_muscle_ache |
      symptoms_new_cough | symptoms_new_sore_throat | symptoms_new_short_breath | symptoms_new_loss_smell)
  ) |>
  # create flag for any symptoms
  dplyr::mutate(
    symptoms_any = (symptoms_fever | symptoms_tiredness | symptoms_headache | symptoms_muscle_ache |
      symptoms_cough | symptoms_sore_throat | symptoms_short_breath | symptoms_loss_smell)
  ) |>
  # create flag for any symptoms
  dplyr::mutate(
    symptoms_count_new_ili = as.integer(symptoms_new_fever) + as.integer(symptoms_new_tiredness) +
      as.integer(symptoms_new_headache) + as.integer(symptoms_new_muscle_ache) + as.integer(symptoms_new_cough) +
      as.integer(symptoms_new_sore_throat) + as.integer(symptoms_new_short_breath) + as.integer(symptoms_new_loss_smell)
  ) |>
  dplyr::mutate(
    window_submission = as.integer(window_submission),
    symptoms_submission = as.integer(symptoms_submission),
    .row = dplyr::row_number()
  ) |>
  age_groups_mutate() |>
  dplyr::mutate(age_group = as.factor(age_group)) |>
  dplyr::filter(rgn21nm != "Scotland") |>
  dplyr::filter(
    submitted_at >= as.Date(study_start_date),
    submitted_at < as.Date(study_end_date)
  )


message("finished load_data script")
