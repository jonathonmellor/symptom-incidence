# script to run Richardson–Lucy deconvolution method to generate incidence
# based on https://pnas.org/doi/full/10.1073/pnas.0902958106

# depends on:
# - prevalence.R
# - duration.R



pacman::p_load(tidyverse, glue, here, dplyr, scales, fastbeta)



ggplot2::theme_set(ggplot2::theme_bw())
source(here::here("utils.R"))


data_output <- fs::dir_create(here::here(
  "outputs",
  "data",
  "sensitivity",
  "incidence",
  "convolution"
))
plot_output <- fs::dir_create(here::here(
  "outputs",
  "figures",
  "sensitivity",
  "incidence",
  "convolution"
))

symptom_opts <- c(
  "cough",
  "fever",
  "sore_throat",
  "short_breath",
  "loss_smell"
)

symptom_labels <- c(
  "cough" = "Cough",
  "fever" = "Fever",
  "sore_throat" = "Sore throat",
  "short_breath" = "Shortness of breath",
  "loss_smell" = "Loss of smell"
)




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
  dplyr::rename(date = submitted_at) |>
  dplyr::arrange(date)

study_start_date <- prevalence_age$date |> min()

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

# a maximum size for the deconvolution
distribution_size <- 60

# generate the draws for the delay convolution for each age group and symptom
duration_age_dist <- duration_age |>
  tidyr::expand_grid(symp_day = seq(0, distribution_size)) |>
  dplyr::mutate(
    delay_pmf = dlnorm(x = symp_day, meanlog = mu, sdlog = sigma),
    .by = c("age_group", ".draw", "symptom")
  )

# the deconvolution is an iterative algorithm, how close to fit it?
tol <- 1e-4


# initialise data frame to fill with each iteration.
incidence_raw <- data.frame()

# extremely slow implementation with loops, but quick to write...
# Alternative would be to nest some dataframes.
# I've done this because it's a one-off sensitivity, and the
# fastbeta::deconvolve doesn't work well with dataframes

# loop over symptom, age group and draw to build the incidence table
for (symptom_i in symptom_opts) {
  message(symptom_i)
  for (age_group_i in unique(duration_age$age_group)) {
    message(age_group_i)
    for (draw_i in seq(min(duration_age_dist$.draw), max(duration_age_dist$.draw))) {
      if (draw_i %% 100 == 0) message(draw_i)

      # select a specific iteration of the duration distribution
      # as per usual, we are matching .draws across data sets,
      # which is arbitrary as they are not jointly modelled.
      dist <- duration_age_dist |>
        dplyr::filter(
          symptom == symptom_i,
          age_group == age_group_i,
          .draw == draw_i
        )

      # take prevalence of a given age/symptom and deconvolve with the
      # duration distribution draw.
      incidence_i <- prevalence_age |>
        dplyr::filter(
          symptom == symptom_i,
          age_group == age_group_i,
          .draw == draw_i
        ) |>
        dplyr::summarise(
          incidence = fastbeta::deconvolve(
            x = prevalence,
            # all symptoms are observed
            prob = 1,
            # pass in the time delay
            delay = dist$delay_pmf,
            # arbitrary tolerate level
            tol = tol
          )$value |>
            # the deconvolution pads the end of the time series, which we don't want
            utils::tail(-distribution_size),
          # keep in population for ease of weighting
          population = unique(population),
          .by = c(".draw", "symptom", "age_group")
        ) |>
        dplyr::mutate(
          t = dplyr::row_number(),
          # we have moved backward in time, but not yet scaled
          incidence = incidence / dist$exp_dur[1]
        )

      # add the data back together
      incidence_raw <- dplyr::bind_rows(
        incidence_raw,
        incidence_i
      )
    }
  }
}

# clean up the data
incidence_regional <- incidence_raw |>
  dplyr::mutate(date = as.Date(study_start_date) + t) |>
  dplyr::select(-t)

# aggregate up with a population weighted sum
incidence_national <- incidence_regional |>
  dplyr::summarise(
    incidence = sum(incidence * population / sum(population)),
    .by = c("date", "symptom", ".draw")
  ) |>
  dplyr::mutate(age_group = "All")

# bring together matching formatting of original incidence
incidence <- dplyr::bind_rows(
  incidence_regional,
  incidence_national
)


# write out
readr::write_csv(incidence,
  file = here::here(
    data_output,
    "symptom_incidence.csv"
  )
)

# read in and summarise
incidence_new_method <- vroom::vroom(
  fs::path(data_output, "symptom_incidence.csv")
) |>
  dplyr::mutate(
    incidence = 1000 * incidence
  ) |>
  dplyr::summarise(
    pi_50 = quantile(incidence, 0.5, na.rm = TRUE),
    pi_95 = quantile(incidence, 0.95, na.rm = TRUE),
    pi_5 = quantile(incidence, 0.05, na.rm = TRUE),
    pi_75 = quantile(incidence, 0.75, na.rm = TRUE),
    pi_25 = quantile(incidence, 0.25, na.rm = TRUE),
    .by = c("symptom", "age_group", "date")
  ) |>
  dplyr::mutate(data = "Deconvolution")


# Process the original WCIS incidence method
# taken from compare.R

# load WCIS incidence
wcis_data <- readr::read_csv(
  file = here::here(
    "outputs",
    "data",
    "symptom_incidence.csv"
  )
) |>
  dplyr::summarise(
    incidence = sum(incidence),
    population = sum(population),
    .by = c("symptom", "date", "age_group", ".draw")
  )

# combine all ages and national into one data frame
wcis_samples <- wcis_data |>
  dplyr::summarise(
    incidence = sum(incidence),
    population = sum(population),
    .by = c("symptom", "date", ".draw")
  ) |>
  dplyr::mutate(age_group = "All") |>
  dplyr::bind_rows(wcis_data)

# calculate quantiles
wcis_summary <- wcis_samples |>
  # convert to a per capita rate
  dplyr::mutate(incidence = 1000 * incidence / population) |>
  dplyr::summarise(
    pi_50 = quantile(incidence, 0.5),
    pi_95 = quantile(incidence, 0.95),
    pi_5 = quantile(incidence, 0.05),
    pi_75 = quantile(incidence, 0.75),
    pi_25 = quantile(incidence, 0.25),
    .by = c("symptom", "age_group", "date")
  ) |>
  dplyr::mutate(data = "Fixed back-shift")


# bring original and deconvolution together
combined_incidence <- dplyr::bind_rows(
  wcis_summary,
  incidence_new_method
) |>
  dplyr::group_by(symptom, age_group, data) |>
  dplyr::mutate(max_date = max(date)) |>
  dplyr::ungroup() |>
  dplyr::group_by(symptom, age_group) |>
  dplyr::filter(date <= min(max_date)) |>
  dplyr::ungroup() |>
  dplyr::select(-max_date)

# make plot for supplement
convolution_plot <- combined_incidence |>
  dplyr::mutate(
    data = factor(data, level = c("Fixed back-shift", "Deconvolution")),
    age_group = factor(age_group, levels = c("All", age_group_order))
  ) |>
  ggplot() +
  geom_line(aes(x = date, y = pi_50, group = data, color = data)) +
  geom_ribbon(aes(x = date, ymax = pi_95, ymin = pi_5, group = data, fill = data), alpha = 0.5) +
  labs(x = "Date", y = "Incidence rate of symptom per day") +
  theme(axis.text.x = element_text(size = 8, angle = 90)) +
  scale_x_date(
    breaks = "3 weeks", minor_breaks = "1 week",
    labels = scales::date_format("%b-%d")
  ) +
  scale_y_continuous(labels = scales::comma) +
  ggh4x::facet_grid2(
    age_group ~ symptom,
    scales = "free_y",
    independent = "y",
    labeller = labeller(symptom = as_labeller(symptom_labels))
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "bottom"
  ) +
  scale_fill_viridis_d(direction = 1, begin = 0.0, end = 0.7) +
  scale_colour_viridis_d(direction = 1, begin = 0.0, end = 0.7)

convolution_plot


# save figures
ggplot2::ggsave(
  here::here(
    plot_output,
    "incidence_plot_deconvolution.png"
  ),
  convolution_plot,
  height = 14,
  width = 12
)

# compare relative difference between the two methods
relative_difference_plot <- combined_incidence |>
  dplyr::mutate(
    data = dplyr::if_else(data == "Deconvolution", "d", "f"),
    age_group = factor(age_group, levels = c("All", age_group_order))
  ) |>
  tidyr::pivot_wider(names_from = data, values_from = dplyr::starts_with("pi")) |>
  dplyr::mutate(
    diff_50 = (pi_50_f - pi_50_d) / pi_50_d,
    diff_5 = (pi_5_f - pi_5_d) / pi_5_d,
    diff_95 = (pi_95_f - pi_95_d) / pi_95_d
  ) |>
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  geom_line(aes(x = date, y = diff_50)) +
  geom_ribbon(aes(x = date, ymax = diff_95, ymin = diff_5), fill = "forestgreen", alpha = 0.5) +
  labs(x = "Date", y = "Relative difference between fixed back-shift and deconvolution (%)") +
  theme(axis.text.x = element_text(size = 8, angle = 90)) +
  scale_x_date(
    breaks = "3 weeks", minor_breaks = "1 week",
    labels = scales::date_format("%b-%d")
  ) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  ggh4x::facet_grid2(
    age_group ~ symptom,
    labeller = labeller(symptom = as_labeller(symptom_labels))
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "bottom"
  )

relative_difference_plot

ggplot2::ggsave(
  here::here(
    plot_output,
    "rel_diff_incidence_plot_deconvolution.png"
  ),
  relative_difference_plot,
  height = 14,
  width = 12
)
