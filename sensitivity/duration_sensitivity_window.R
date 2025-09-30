# Estimate symptom duration distribution using Epidist by age group.
# Sensitivity for how wide the duration window is on the secondary event.

# depends on:
# - load_data.R

# install epidist if not already installed
if (!require("epidist")) {
  remotes::install_github(
    file.path("epinowcast", "epidist"),
    dependencies = TRUE
  )
}

pacman::p_load(tidyverse, glue, epidist, dplyr, modelr)

set.seed(13)

source(here::here("syndromic_surveillance", "manuscripts", "symptom_incidence", "load_data.R"))

fs::dir_create(here::here(
  "outputs",
  "figures",
  "sensitivity",
  "duration",
  "window")
)

fs::dir_create(here::here(
  "outputs",
  "data",
  "sensitivity",
  "duration",
  "window")
)

ggplot2::theme_set(ggplot2::theme_bw())

# selected symptoms to calculate durations for
symptom_opts <- c(
  "cough",
  "fever",
  "sore_throat",
  "short_breath",
  "loss_smell"
)

# chosen a range of values on the scale of ILI infections (1-2 weeks).
# This parameter is an addition to the window length of the upper bound, so should be
# a smaller scale than the total infection. Have selected a coarse grid for ease
# of interpretation
window_extension_opts <- c(
  0, 3, 7,
  10)

symptom_labels <- c(
  "cough" = "Cough",
  "fever" = "Fever",
  "sore_throat" = "Sore throat",
  "short_breath" = "Shortness of breath",
  "loss_smell" = "Loss of smell"
)


# data preprocessing
symptoms_data <- symptoms |>
  # rename to symptoms_ongoing for pivot longer
  dplyr::rename_with(
    .fn = \(col_name) stringr::str_replace(col_name, "symptoms_", "symptoms_ongoing_"),
    .cols = dplyr::starts_with("symptoms_") & !dplyr::starts_with("symptoms_new") & dplyr::contains(symptom_opts)
  ) |>
  dplyr::select(
    submitted_at,
    participant_id,
    symptoms_onset,
    age_group,
    starts_with("symptoms_new"),
    starts_with("symptoms_ongoing")
  ) |>
  dplyr::rename_with(~ stringr::str_remove(., "symptoms_")) |>
  tidyr::pivot_longer(
    cols = c(dplyr::starts_with("new"), dplyr::starts_with("ongoing")),
    names_to = c(".value", "symptom"),
    names_pattern = "(.*?)_(.*)"
  ) |>
  # filter to only selected symptoms
  dplyr::filter(stringr::str_detect(symptom, paste(symptom_opts, collapse = "|"))) |>
  # filter to only "new" symptoms
  dplyr::filter(new == TRUE) |>
  # truncate durations longer than 60 days
  dplyr::filter(submitted_at + 1 - onset < 60)


# raw data for comparison
raw_data <- symptoms_data |>
  dplyr::mutate(duration = submitted_at + 1 - onset) |>
  dplyr::group_by(symptom, age_group, duration) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::ungroup()

for (window_extension in window_extension_opts) {
  for (symptom_name in symptom_opts) {
    message(glue::glue("running {symptom_name}, window: {window_extension}"))

    # get data in correct form
    data_prep <- symptoms_data |>
      dplyr::filter(symptom == symptom_name) |>
      # convert data to correct format for epidist
      dplyr::mutate(
        pdate_lwr = onset,
        pdate_upr = pdate_lwr + 1,
        sdate_lwr = submitted_at + 1,
        sdate_upr = sdate_lwr + 1 + window_extension
      ) |>
      dplyr::summarise(n = dplyr::n(),
        .by = c("pdate_lwr", "pdate_upr",
          "sdate_lwr", "sdate_upr", "age_group")) |>
      epidist::as_epidist_aggregate_data() |>
      epidist::as_epidist_marginal_model()

    # run model with age group as a covariate
    fit <- data_prep |>
      epidist::epidist(
        formula = brms::bf(mu ~ 1 + age_group,
          sigma ~ 1 + age_group),
        algorithm = "sampling",
        iter = 1000,
        cores = 2,
        chains = 2,
        backend = "cmdstanr"
      )

    # generate samples
    draws <- data_prep |>
      modelr::data_grid(age_group) |>
      dplyr::mutate(relative_obs_time = NA, pwindow = NA, swindow = NA,
        delay_upr = NA) |>
      tidybayes::add_epred_draws(fit, dpar = TRUE) |>
      dplyr::ungroup() |>
      dplyr::select(age_group, .draw, mu, sigma) |>
      dplyr::mutate(symptom = symptom_name,
        window_extension = window_extension)

    # mean duration by age group
    print(draws |>
      dplyr::mutate(exp_dur = exp(mu + sigma^2 / 2)) |>
      dplyr::summarise(mean_dur = mean(exp_dur), .by = age_group))


    # write up
    readr::write_csv(draws,
      file = here::here(
        "outputs",
        "data",
        "sensitivity",
        "duration",
        "window",
        glue::glue("symptom_durations_samples_{symptom_name}_{window_extension}.csv")
      )
    )
  }


  # combine outputs into a single data frame
  draws_cmb <- purrr::imap(as.list(symptom_opts), function(symptom_name, ind) {
    readr::read_csv(
      file = here::here(
        "outputs",
        "data",
        "sensitivity",
        "duration",
        "window",
        glue::glue("symptom_durations_samples_{symptom_name}_{window_extension}.csv")
      )
    )
  }) |>
    dplyr::bind_rows()

  symptom_durations_cmb <- tidyr::expand_grid(
    draws_cmb,
    duration = seq(0, 100, by = 0.5)
  ) |>
    dplyr::mutate(density = dlnorm(duration, mu, sigma))

  # generate figure by age group
  dur_plot_cmb <- symptom_durations_cmb |>
    dplyr::summarise(
      pi_50 = quantile(density, 0.5),
      pi_95 = quantile(density, 0.95),
      pi_5 = quantile(density, 0.05),
      pi_75 = quantile(density, 0.75),
      pi_25 = quantile(density, 0.25),
      .by = c("symptom", "age_group", "duration")
    ) |>
    ggplot() +
    geom_ribbon(aes(x = duration, ymax = pi_95, ymin = pi_5, fill = symptom, group = symptom),
      alpha = 0.5) +
    geom_line(aes(x = duration, y = pi_50, color = symptom, group = symptom)) +
    xlim(0, 100) +
    labs(x = "Duration (days)", y = "Probability density") +
    xlim(0, 30) +
    facet_grid(
      factor(age_group, age_group_order) ~ symptom,
      scales = "fixed",
      labeller = labeller(symptom = as_labeller(symptom_labels))
    ) +
    theme(legend.position = "none")

  # save figure
  ggplot2::ggsave(here::here(
    "outputs",
    "figures",
    "sensitivity",
    "duration",
    "window",
    glue::glue("dur_plot_cmb_{window_extension}.png")
  ),
  dur_plot_cmb,
  height = 10,
  width = 8)

  age_group_labels <- c("3-17", "18-34", "35-44", "45-54", "55-64", "65-74", "75+")

  # compare mean durations across symptoms and age groups
  mean_dur_plot <- draws_cmb |>
    dplyr::mutate(exp_dur = exp(mu + sigma^2 / 2)) |>
    dplyr::summarise(
      .mean = mean(exp_dur),
      .q50 = quantile(exp_dur, 0.5),
      .q95 = quantile(exp_dur, 0.95),
      .q5 = quantile(exp_dur, 0.05),
      .q75 = quantile(exp_dur, 0.75),
      .q25 = quantile(exp_dur, 0.25), .by = c("age_group", "symptom")
    ) |>
    ggplot(aes(x = factor(age_group, age_group_order), color = symptom, group = symptom)) +
    geom_linerange(aes(ymin = .q5, ymax = .q95)) +
    geom_linerange(aes(ymin = .q25, ymax = .q75), linewidth = 1) +
    geom_point(aes(y = .q50)) +
    scale_x_discrete(labels = age_group_labels) +
    theme(axis.text.x = element_text(size = 8)) +
    labs(y = "Mean duration (days)", x = "Age group") +
    facet_grid(~symptom, labeller = as_labeller(symptom_labels)) +
    theme(legend.position = "none")

  # save figure
  ggplot2::ggsave(here::here(
    "outputs",
    "figures",
    "sensitivity",
    "duration",
    "window",
    glue::glue("dur_plot_mean_{window_extension}.png")
  ),
  mean_dur_plot,
  height = 3,
  width = 14)

}
message("DURATION MODELLING DONE")
