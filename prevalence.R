# script to model WCIS symptom prevalence.
pacman::p_load(tidyverse, mgcv, gratia, glue)

# depends on:
# - load_data.R


source(here::here("utils.R"))
source(here::here("load_data.R"))

fs::dir_create(here::here("outputs", "data"))

study_start_date <- "2023-11-14"
study_end_date <- "2024-03-07"


outcomes <- c(
  "symptoms_new_cough",
  "symptoms_new_fever",
  "symptoms_new_sore_throat",
  "symptoms_new_short_breath",
  "symptoms_new_loss_smell"
)


for (outcome in outcomes) {
  message(glue::glue("running: {outcome}"))

  symptoms_model_ready <- symptoms |>
    dplyr::filter(rgn21nm != "Scotland") |>
    dplyr::mutate(
      age_group = factor(age_group, levels = age_group_order),
      region_name = factor(rgn21nm, levels = region_order)
    ) |>
    dplyr::filter(submitted_at >= study_start_date & submitted_at < as.Date(study_end_date)) |>
    dplyr::select(all_of(c("window_submission", "submitted_at", "age_group", "region_name", outcome))) |>
    # define model time indexed on submission date
    dplyr::mutate(t = 1 + as.integer(submitted_at) - as.integer(min(submitted_at))) |>
    dplyr::summarise(
      n_trials = dplyr::n(),
      n_outcome = sum(!!sym(outcome)),
      .by = c("age_group", "region_name", "window_submission", "t", "submitted_at")
    )


  # calculate how present each testing day is to create a weighting on the modelled
  # effects
  window_weighting <- symptoms_model_ready |>
    dplyr::summarize(window_n = n(), .by = "window_submission") |>
    dplyr::mutate(window_weight = window_n / sum(window_n)) |>
    dplyr::select(all_of(c("window_submission", "window_weight")))

  # create set of unique combinations of the modelled data
  # to perform sample generation and postrat on
  symptoms_unique <- tidyr::expand_grid(
    submitted_at = symptoms_model_ready$submitted_at |> unique(),
    age_group = symptoms_model_ready$age_group |> unique(),
    window_submission = symptoms_model_ready$window_submission |> unique(),
    region_name = symptoms_model_ready$region_name |> unique()
  ) |>
    # define model time indexed on submission date
    dplyr::mutate(t = 1 + as.integer(submitted_at) - as.integer(min(submitted_at))) |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    # need at prediction time, but just a column rather than a value
    dplyr::mutate(n_trials = 0)

  # table of true population used to stratify on
  poststrat_age <- poststrat |>
    age_groups_mutate() |>
    dplyr::mutate(age_group = factor(age_group, levels = age_group_order)) |>
    dplyr::summarise(population = sum(population), .by = c("age_group", "region_name"))


  message("model fitting")


  # prevalence over time, indexed by the submission date
  model <- mgcv::gam(
    formula = glue::glue("cbind(n_outcome,n_trials) ~
  window_submission +
  s(t, bs='gp', m=c(-2, 50, 2)) +
  s(t, by=age_group, bs='gp', m=c(-2, 20, 2)) +
  s(region_name, bs='re')") |> as.formula(),
    data = symptoms_model_ready,
    family = "binomial(logit)"
  )


  message("post processing")


  samples <- gratia::fitted_samples(model,
    data = symptoms_unique,
    scale = "linear_predictor",
    method = "mh",
    seed = 1,
    n = 1000
  ) |>
    # convert from linear predictor scale to a probability
    dplyr::mutate(.p = 1 / (1 + exp(-.fitted))) |>
    dplyr::select(-.fitted) |>
    # bring back covariates into the samples
    dplyr::right_join(symptoms_unique, by = ".row") |>
    dplyr::mutate(age_group = factor(age_group, levels = age_group_order)) |>
    # reweights over the window effect by how much each window occurs
    dplyr::left_join(window_weighting, by = "window_submission") |>
    dplyr::summarise(.p = sum(.p * window_weight), .by = c("submitted_at", ".draw", "region_name", "age_group")) |>
    # bring in postratification information
    dplyr::left_join(poststrat_age, by = c("age_group", "region_name")) |>
    dplyr::mutate(.p_n = .p * population)

  # MRP age prevalence
  predictions_age <- samples |>
    dplyr::summarise(
      .p = sum(.p_n) / sum(population),
      .by = c("submitted_at", ".draw", "age_group")
    ) |>
    dplyr::summarise(
      pi_50 = quantile(.p, 0.5),
      pi_95 = quantile(.p, 0.95),
      pi_5 = quantile(.p, 0.05),
      pi_75 = quantile(.p, 0.75),
      pi_25 = quantile(.p, 0.25),
      .by = c("submitted_at", "age_group")
    ) |>
    dplyr::mutate(
      outcome = outcome,
      breakdown = "age"
    )

  # MRP regional prevalence
  predictions_region <- samples |>
    dplyr::summarise(
      .p = sum(.p_n) / sum(population),
      .by = c("submitted_at", ".draw", "region_name")
    ) |>
    dplyr::summarise(
      pi_50 = quantile(.p, 0.5),
      pi_95 = quantile(.p, 0.95),
      pi_5 = quantile(.p, 0.05),
      pi_75 = quantile(.p, 0.75),
      pi_25 = quantile(.p, 0.25),
      .by = c("submitted_at", "region_name")
    ) |>
    dplyr::mutate(
      outcome = outcome,
      breakdown = "region"
    )

  # MRP national prevalence, weighted
  predictions_national <- samples |>
    dplyr::summarise(.p = sum(.p_n) / sum(population), .by = c("submitted_at", ".draw")) |>
    dplyr::summarise(
      pi_50 = quantile(.p, 0.5),
      pi_95 = quantile(.p, 0.95),
      pi_5 = quantile(.p, 0.05),
      pi_75 = quantile(.p, 0.75),
      pi_25 = quantile(.p, 0.25),
      .by = c("submitted_at")
    ) |>
    dplyr::mutate(
      outcome = outcome,
      breakdown = "nation"
    )

  predictions <- dplyr::bind_rows(
    predictions_age,
    predictions_region,
    predictions_national
  )

  # WRITE OUT ####
  # summary stats
  readr::write_csv(predictions,
    file = here::here(
      "outputs",
      "data",
      glue::glue("prevalence_summary_{outcome}.csv")
    )
  )
  # samples
  readr::write_csv(samples,
    file = here::here(
      "outputs",
      "data",
      glue::glue("prevalence_samples_{outcome}.csv")
    )
  )

  message(glue::glue("done writing prevalence {outcome}"))
}

message("PREVALENCE MODELLING DONE")
