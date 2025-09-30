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

library(ggplot2)
ggplot2::theme_set(ggplot2::theme_bw())

source(here::here("utils.R"))

formulas <- list(
  "mu" = brms::bf(mu ~ 1 + age_group),
  "mu_sigma" = brms::bf(
    mu ~ 1 + age_group,
    sigma ~ 1 + age_group
  )
)

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
    "formula_name" = names(formulas)
  ) |> as.list(),
  function(symptom_name, formula_name) {
    readr::read_csv(
      file = here::here(
        "outputs",
        "data",
        "sensitivity",
        "duration",
        "sigma_analysis",
        glue::glue("symptom_durations_samples_{symptom_name}_{formula_name}.csv")
      )
    ) |> dplyr::mutate(formula_name = formula_name)
  }
) |>
  dplyr::bind_rows()

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
    .q25 = quantile(exp_dur, 0.25), .by = c("age_group", "symptom", "formula_name")
  ) |>
  dplyr::mutate(formula_name = factor(formula_name)) |>
  ggplot(aes(x = factor(age_group, age_group_order), color = formula_name, group = formula_name)) +
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
