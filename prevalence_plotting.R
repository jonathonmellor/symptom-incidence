# plot the prevalence results used for the supplimentary material
pacman::p_load(tidyverse, mgcv, glue, gratia, patchwork, forcats)

source(here::here("utils.R"))

fs::dir_create(here::here(
  "outputs",
  "figures",
  "prevalence"
))

ggplot2::theme_set(ggplot2::theme_bw())

outcomes <- c(
  "symptoms_new_cough",
  "symptoms_new_fever",
  "symptoms_new_sore_throat",
  "symptoms_new_short_breath",
  "symptoms_new_loss_smell"
)

symptom_labels <- c(
  "cough" = "Cough",
  "fever" = "Fever",
  "sore_throat" = "Sore throat",
  "short_breath" = "Shortness of breath",
  "loss_smell" = "Loss of smell"
)

all_outcomes <- data.frame()

for (outcome in outcomes) {
  predictions <- readr::read_csv(
    file = here::here(
      "outputs",
      "data",
      glue::glue("prevalence_summary_{outcome}.csv")
    )
  )
  all_outcomes <- dplyr::bind_rows(all_outcomes, predictions)
}

prevalence <- all_outcomes |>
  dplyr::mutate(symptom = stringr::str_remove(outcome, "symptoms_new_"))

age_plot <- prevalence |>
  dplyr::filter(breakdown %in% c("age", "nation")) |>
  dplyr::mutate(age_group = ifelse(breakdown == "nation", "All", age_group)) |>
  ggplot() +
  geom_ribbon(aes(
    x = submitted_at, ymax = 100 * pi_95, ymin = 100 * pi_5,
    fill = symptom, group = symptom
  ), alpha = 0.6) +
  geom_line(aes(x = submitted_at, y = 100 * pi_50, group = symptom)) +
  ylab("Prevalence of symptom (%)") +
  xlab("Submitted date") +
  theme(axis.text.x = element_text(size = 8, angle = 90)) +
  scale_x_date(
    breaks = "2 weeks", minor_breaks = "4 week",
    labels = scales::date_format("%b-%d")
  ) +
  ggh4x::facet_grid2(
    factor(age_group, c("All", age_group_order)) ~ symptom,
    scales = "free_y",
    independent = "y",
    labeller = labeller(symptom = as_labeller(symptom_labels))
  ) +
  theme(legend.position = "none")

age_plot

ggplot2::ggsave(
  here::here(
    "outputs", "figures",
    "prevalence",
    "age_prev.png"
  ),
  age_plot,
  height = 12,
  width = 12
)

region_plot <- prevalence |>
  dplyr::filter(breakdown == "region") |>
  ggplot() +
  geom_ribbon(aes(
    x = submitted_at, ymax = 100 * pi_95, ymin = 100 * pi_5,
    fill = symptom, group = symptom
  ), alpha = 0.6) +
  geom_line(aes(x = submitted_at, y = 100 * pi_50, group = symptom)) +
  ylab("proportion of population with symptoms (%)") +
  xlab("submitted date") +
  theme(axis.text.x = element_text(size = 5, angle = 90)) +
  scale_x_date(
    breaks = "2 weeks", minor_breaks = "4 week",
    labels = scales::date_format("%b-%d")
  ) +
  facet_wrap(~region_name, scales = "fixed") +
  theme(legend.position = "bottom")

region_plot

ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "prevalence",
    "region_prev.png"
  ),
  region_plot,
  height = 8,
  width = 12
)

nation_plot <- prevalence |>
  dplyr::filter(breakdown == "nation") |>
  ggplot() +
  geom_ribbon(aes(
    x = submitted_at, ymax = 100 * pi_95, ymin = 100 * pi_5,
    fill = symptom, group = symptom
  ), alpha = 0.5) +
  geom_line(aes(x = submitted_at, y = 100 * pi_50, group = symptom)) +
  ylab(glue::glue("Proportion of population with symptoms (%)")) +
  xlab("Submitted date") +
  facet_wrap(~symptom, scales = "free_y", , labeller = labeller(symptom = as_labeller(symptom_labels))) +
  scale_x_date(breaks = "2 weeks", minor_breaks = "1 week", labels = scales::date_format("%b-%d")) +
  theme(legend.position = "bottom")

nation_plot

ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "prevalence",
    "nation_prev.png"
  ),
  nation_plot,
  height = 8,
  width = 12
)
