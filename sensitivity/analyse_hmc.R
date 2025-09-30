# script to compare the results of HMC and the different approximation approaches
# for `epidist`. Decision: go with HMC using the marginal model.

# relies on `duration_sensitivity_hmc.R` being run.


library(ggplot2)

# selected symptoms to calculate durations for
symptom_opts <- c(
  "cough",
  "fever",
  "sore_throat",
  "short_breath",
  "loss_smell"
)

source(here::here("utils.R"))

age_group_labels <- c("3-17", "18-34", "35-44", "45-54", "55-64", "65-74", "75+")

symptom_labels <- c(
  "cough" = "Cough",
  "fever" = "Fever",
  "sore_throat" = "Sore throat",
  "short_breath" = "Shortness of breath",
  "loss_smell" = "Loss of smell"
)

draws_cmb_hmc <- purrr::pmap(tidyr::expand_grid("symptom_name" = symptom_opts) |> as.list(),
  function(symptom_name) {
    readr::read_csv(
      file = here::here(
        "outputs",
        "data",
        "sensitivity",
        "duration",
        "hmc",
        glue::glue("symptom_durations_samples_{symptom_name}.csv")
      )
    ) |> dplyr::mutate(sampling_method = "sigma varied")
  }) |>
  dplyr::bind_rows()

draws_cmb_approx <- purrr::pmap(
  tidyr::expand_grid("symptom_name" = symptom_opts) |> as.list(),
  function(symptom_name) {
    readr::read_csv(
      file = here::here(
        "outputs",
        "data",
        glue::glue("symptom_durations_samples_{symptom_name}.csv")
      )
    ) |> dplyr::mutate(sampling_method = "sigma fixed")
  }) |>
  dplyr::bind_rows()

draws_cmb <- dplyr::bind_rows(
  draws_cmb_hmc,
  draws_cmb_approx
)

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
    .q25 = quantile(exp_dur, 0.25), .by = c("age_group", "symptom", "sampling_method")
  ) |>
  ggplot(aes(x = factor(age_group, age_group_order), color = sampling_method, group = sampling_method)) +
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

ggplot2::ggsave(here::here(
  "outputs",
  "figures",
  "sensitivity",
  "duration",
  "window",
  glue::glue("dur_plot_mean_all_windows.png")
),
mean_dur_plot,
height = 8,
width = 14)



# generate figure by age group
dur_plot_cmb <- symptom_durations_cmb |>
  dplyr::summarise(
    pi_50 = quantile(density, 0.5),
    pi_95 = quantile(density, 0.95),
    pi_5 = quantile(density, 0.05),
    pi_75 = quantile(density, 0.75),
    pi_25 = quantile(density, 0.25),
    .by = c("symptom", "age_group", "duration", "sampling_method")
  ) |>
  ggplot() +
  geom_ribbon(aes(x = duration, ymax = pi_95, ymin = pi_5, fill = sampling_method, group = sampling_method),
    alpha = 0.5) +
  geom_line(aes(x = duration, y = pi_50, color = sampling_method, group = sampling_method)) +
  xlim(0, 100) +
  labs(x = "Duration (days)", y = "Probability density") +
  xlim(0, 30) +
  facet_grid(
    factor(age_group, age_group_order) ~ symptom,
    scales = "fixed",
    labeller = labeller(symptom = as_labeller(symptom_labels))
  ) +
  theme(legend.position = "none")

dur_plot_cmb
