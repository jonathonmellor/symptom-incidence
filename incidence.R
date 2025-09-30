# script to model symptom incidence from WCIS data

# depends on:
# - prevalence.R
# - duration.R


pacman::p_load(tidyverse, glue, here, dplyr, scales)

# note: takes ~8 mins
source(here::here("prevalence.R"))
# note: takes ~1.5 hours
source(here::here("duration.R"))

source(here::here("utils.R"))

ggplot2::theme_set(ggplot2::theme_bw())

fs::dir_create(here::here(
  "outputs",
  "figures",
  "incidence"
))

symptom_opts <- c(
  "cough",
  "fever",
  "sore_throat",
  "short_breath",
  "loss_smell"
)


study_start_date <- "2023-11-14"
study_end_date <- "2024-03-07"

# get prevalence
prevalence_age <- purrr::imap(as.list(symptom_opts), function(symptom_name, ind) {
  readr::read_csv(
    file = here::here(
      "outputs",
      "data",
      glue::glue("prevalence_samples_symptoms_new_{symptom_name}.csv")
    ),
    show_col_types = FALSE
  ) |>
    dplyr::mutate(symptom = as.factor(symptom_name))
}) |>
  # combine into a single data frame
  dplyr::bind_rows() |>
  dplyr::summarise(
    prevalence = sum(.p_n) / sum(population),
    population = sum(population),
    .by = c("symptom", "submitted_at", ".draw", "age_group")
  ) |>
  dplyr::rename(date = submitted_at)

# get durations
duration_age <- purrr::imap(as.list(symptom_opts), function(symptom_name, ind) {
  readr::read_csv(
    file = here::here(
      "outputs",
      "data",
      glue::glue("symptom_durations_samples_{symptom_name}.csv")
    )
  ) |>
    dplyr::mutate(symptom = symptom_name)
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(exp_dur = exp(mu + sigma^2 / 2))

# calculate incidence for each draw
incidence <- dplyr::left_join(prevalence_age,
  duration_age,
  by = c("symptom", "age_group", ".draw")
) |>
  dplyr::mutate(
    incidence = prevalence * population / exp_dur,
    # shift incidence back in time, since prevalence is by submission date
    min_date = min(date) - round(min(exp_dur) / 2),
    max_date = max(date) - round(max(exp_dur) / 2),
    date = date - round(exp_dur / 2),
    .by = c("symptom")
  ) |>
  dplyr::filter((date >= min_date) & (date <= max_date))

# calculate national incidence
incidence_national <- incidence |>
  dplyr::summarise(
    incidence = sum(incidence),
    population = sum(population),
    .by = c("symptom", "date", ".draw")
  )

# plot incidence by age
incidence_plot <- incidence |>
  dplyr::summarise(
    incidence_q50 = quantile(incidence, 0.5),
    incidence_q95 = quantile(incidence, 0.95),
    incidence_q5 = quantile(incidence, 0.05),
    incidence_q75 = quantile(incidence, 0.75),
    incidence_q25 = quantile(incidence, 0.25),
    .by = c("symptom", "date", "age_group")
  ) |>
  dplyr::mutate(symptom = sub("_", " ", symptom)) |>
  ggplot() +
  geom_ribbon(aes(x = date, ymax = incidence_q95, ymin = incidence_q5, fill = symptom, group = symptom),
    alpha = 0.6
  ) +
  geom_line(aes(x = date, y = incidence_q50, group = symptom)) +
  ylab("incidence") +
  xlab("date") +
  theme(axis.text.x = element_text(size = 8, angle = 90)) +
  scale_x_date(
    breaks = "2 weeks", minor_breaks = "4 week",
    labels = scales::date_format("%b-%d")
  ) +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(symptom ~ factor(age_group, age_group_order), scales = "free_y") +
  theme(legend.position = "bottom")

# plot national incidence
incidence_national_plot <- incidence_national |>
  dplyr::summarise(
    incidence_q50 = quantile(incidence, 0.5),
    incidence_q95 = quantile(incidence, 0.95),
    incidence_q5 = quantile(incidence, 0.05),
    incidence_q75 = quantile(incidence, 0.75),
    incidence_q25 = quantile(incidence, 0.25),
    .by = c("symptom", "date")
  ) |>
  dplyr::mutate(symptom = sub("_", " ", symptom)) |>
  ggplot() +
  geom_ribbon(aes(x = date, ymax = incidence_q95, ymin = incidence_q5, fill = symptom, group = symptom),
    alpha = 0.6
  ) +
  geom_line(aes(x = date, y = incidence_q50, group = symptom)) +
  ylab("incidence") +
  xlab("date") +
  theme(axis.text.x = element_text(size = 8, angle = 90)) +
  scale_x_date(
    breaks = "2 weeks", minor_breaks = "4 week",
    labels = scales::date_format("%b-%d")
  ) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~symptom, scales = "free_y") +
  theme(legend.position = "bottom")

# write out
readr::write_csv(incidence,
  file = here::here(
    "outputs",
    "data",
    "symptom_incidence.csv"
  )
)

# save figures
ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "incidence",
    "incidence_plot_age.png"
  ),
  incidence_plot,
  height = 8,
  width = 12
)

ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "incidence",
    "incidence_plot_nation.png"
  ),
  incidence_national_plot,
  height = 8,
  width = 12
)

message("INCIDENCE MODELLING DONE")
