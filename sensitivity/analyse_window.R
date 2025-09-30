# script to explore different window lengths of upper bound duration.
# This code is a near direct copy of `duration.R` with additional parameter
# for the window_extension on the secondary event upper limit.

# relies on `duration_sensitivity.R` being run.

# selected symptoms to calculate durations for
symptom_opts <- c(
  "cough",
  "fever",
  "sore_throat",
  "short_breath",
  "loss_smell"
)

source(here::here("utils.R"))

library(ggplot2)
ggplot2::theme_set(ggplot2::theme_bw())


window_extension_opts <- c(0, 3, 7, 10)

age_group_labels <- c("3-17", "18-34", "35-44", "45-54", "55-64", "65-74", "75+")

symptom_labels <- c(
  "cough" = "Cough",
  "fever" = "Fever",
  "sore_throat" = "Sore throat",
  "short_breath" = "Shortness of breath",
  "loss_smell" = "Loss of smell"
)

draws_cmb <- purrr::pmap(
  tidyr::expand_grid(
    "symptom_name" = symptom_opts,
    "window_extension" = window_extension_opts
  ) |> as.list(),
  function(symptom_name, window_extension) {
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
  }
) |>
  dplyr::bind_rows() |>
  # convergence issues with max delay so removed.
  dplyr::filter(window_extension != 10)

symptom_durations_cmb <- tidyr::expand_grid(
  draws_cmb,
  duration = seq(0, 100, by = 0.5)
) |>
  dplyr::mutate(density = dlnorm(duration, mu, sigma))


mean_dur_plot <- draws_cmb |>
  dplyr::mutate(exp_dur = exp(mu + sigma^2 / 2)) |>
  dplyr::summarise(
    .mean = mean(exp_dur),
    .q50 = quantile(exp_dur, 0.5),
    .q95 = quantile(exp_dur, 0.95),
    .q5 = quantile(exp_dur, 0.05),
    .q75 = quantile(exp_dur, 0.75),
    .q25 = quantile(exp_dur, 0.25), .by = c("age_group", "symptom", "window_extension")
  ) |>
  dplyr::mutate(`window extension` = factor(window_extension)) |>
  ggplot(aes(x = factor(age_group, age_group_order), color = `window extension`, group = `window extension`)) +
  geom_linerange(aes(ymin = .q5, ymax = .q95), position = position_dodge2(width = 0.5)) +
  geom_linerange(aes(ymin = .q25, ymax = .q75), linewidth = 1, position = position_dodge2(width = 0.5)) +
  geom_point(aes(y = .q50), position = position_dodge2(width = 0.5)) +
  scale_x_discrete(labels = age_group_labels) +
  theme(axis.text.x = element_text(size = 8)) +
  labs(y = "Mean duration (days)", x = "Age group") +
  facet_grid(~symptom, labeller = as_labeller(symptom_labels)) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 30))

mean_dur_plot

ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "sensitivity",
    "duration",
    "window",
    glue::glue("dur_plot_mean_all_windows.png")
  ),
  mean_dur_plot,
  height = 8,
  width = 14
)
