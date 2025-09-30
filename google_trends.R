# script to model google symptoms data adjusting for DoW effects

# need version >= 0.10.0 for growth rates
stopifnot(packageVersion("gratia") >= "0.10.0")

pacman::p_load(tidyverse, mgcv, glue, gratia, patchwork, forcats)

box::use(
  # DATABASE ACCESS - REMOVED
  prj / projection_plots
)

# script to load all needed data

message("loading data")

db <- # DATABASE ACCESS - REMOVED

study_start_date <- "2023-11-14"
study_end_date <- "2024-03-07"


gt_symptoms <- db$table_name |>
  dplyr::collect() |>
  dplyr::select(date, symptom, value) |>
  dplyr::filter(date >= as.Date(study_start_date),
    date < as.Date(study_end_date))

message("finished loading script")

fs::dir_create(here::here(
  "outputs",
  "data")
)


ggplot2::theme_set(projection_plots$theme_pancasts())

gt_symptom_opts <- c("fever",
  "headache",
  "fatigue",
  "cough",
  "sore_throat",
  "shortness_of_breath",
  "common_cold",
  "myalgia",
  "anosmia"
)

dotw_order <- c(
  "Sunday",
  "Monday",
  "Tuesday",
  "Wednesday",
  "Thursday",
  "Friday",
  "Saturday")


input_data <- gt_symptoms |>
  dplyr::filter(symptom %in% gt_symptom_opts) |>
  dplyr::arrange(symptom, date) |>
  dplyr::mutate(dotw = weekdays(date)) |>
  dplyr::mutate(is_bank_holiday = dplyr::case_when(
    date == "2023-12-25" ~ 1,
    date == "2023-12-26" ~ 1,
    .default = 0
  )) |>
  dplyr::mutate(
    t = 1 + as.integer(date) - as.integer(min(date)),
    symptom = as.factor(symptom),
    dotw = as.integer(factor(dotw, levels = dotw_order)),
    .by = "symptom") |>
  # we are separating by symptom when we model later so want a row per symptom
  dplyr::mutate(.row = dplyr::row_number(), .by = "symptom")

message("Fitting models")

models <- input_data |>
  tidyr::nest(data = -symptom) |>
  # map over each symptom to fit a model independently
  dplyr::mutate(
    model = purrr::map(data,
      \(data) {
        # divide value by 100 to scale to (0,1)
        mgcv::gam(value / 100 ~
          1 +
          # set a lengthscale of 7 as small relative to
          # the timeseries length
          s(t, bs = "gp", m = c(2, 7, 2))  +
          # choosing k=4 as a minimum
          s(dotw, bs = "cc", k = 4) +
          is_bank_holiday,
        data = data,
        method = "REML",
        family = betar(link = "logit"))
      })
  )

message("Generating natural scale samples")


samples <- models |>
  # use the data from the model (not passing in data)
  # as we are working within sample (and makes mapping easier)
  dplyr::mutate(posterior_samples = purrr::map(
    model,
    \(m) {
      gratia::fitted_samples(m,
        scale = "response",
        method = "mh",
        seed = 1,
        # exclude DOW and holiday to show trend rather
        # than artifacts
        exclude = c("s(dotw)", "is_bank_holiday"
        ),
        n = 1000)
    }
  )) |>
  dplyr::select(-c(data, model)) |>
  tidyr::unnest("posterior_samples") |>
  dplyr::right_join(input_data, by = c("symptom", ".row")) |>
  # convert fitted values back to RV scale
  dplyr::mutate(.fitted =  100 * .fitted)

# summarise model predictions
results <- samples |>
  dplyr::summarise(
    pi_50 = quantile(.fitted, 0.5),
    pi_95 = quantile(.fitted, 0.95),
    pi_5 = quantile(.fitted, 0.05),
    pi_75 = quantile(.fitted, 0.75),
    pi_25 = quantile(.fitted, 0.25),
    value = unique(value),
    .by = c("date", "symptom"))

# make a nice plot
natural_plot <- results |>
  dplyr::mutate(symptom = stringr::str_replace_all(symptom, "_", " ")) |>
  ggplot2::ggplot() +
  geom_point(aes(x = date, y = value), alpha = .5, size = 0.5) +
  geom_line(aes(x = date, y = pi_50, color = "central estimate")) +
  geom_ribbon(aes(x = date, ymin = pi_5, ymax = pi_95, fill = "90% interval"), alpha = 0.3) +
  facet_wrap(~symptom, scales = "free_y") +
  scale_color_manual(name = "", values = c("central estimate" = "forestgreen")) +
  scale_fill_manual(name = "", values = c("90% interval" = "forestgreen")) +
  ggtitle("Symptom Trends over Winter 2023/24", subtitle = "Day of week adjusted, gaussian process smoother") +
  ylab("relative search volume")

natural_plot


message("Generating growth-rate scale samples")

gr_samples <- models |>
  dplyr::mutate(gr_posterior_samples = purrr::pmap(
    list(m = model, d = data),
    \(m, d) {
      gratia::derivative_samples(object = m,
        data = d,
        focal = "t",
        type = "central",
        # doing this on the linear predictor scale gives us our `log` scale
        scale = "linear_predictor",
        # we aren't so worried about showing the day of week
        exclude = c("s(dotw)", "is_bank_holiday"),
        n_sim = 1000,
        eps = 0.001)
    }
  )) |>
  dplyr::select(-c(data, model)) |>
  tidyr::unnest("gr_posterior_samples") |>
  dplyr::select(-.focal) |>
  dplyr::right_join(input_data, by = c("symptom", ".row")) |>
  # convert to a %
  dplyr::mutate(.derivative = 100 * .derivative)

gr_results <- gr_samples |>
  dplyr::summarise(
    pi_50 = quantile(.derivative, 0.5),
    pi_95 = quantile(.derivative, 0.95),
    pi_5 = quantile(.derivative, 0.05),
    pi_75 = quantile(.derivative, 0.75),
    pi_25 = quantile(.derivative, 0.25),
    .by = c("date", "symptom"))

gr_plot <- gr_results |>
  dplyr::mutate(symptom = stringr::str_replace_all(symptom, "_", " ")) |>
  ggplot2::ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(x = date, y = pi_50, color = "central estimate")) +
  geom_ribbon(aes(x = date, ymax = pi_95, ymin = pi_5, fill = "90% interval"), alpha = 0.3) +
  facet_wrap(~symptom, scales = "free_y") +
  scale_color_manual(name = "", values = c("central estimate" = "forestgreen")) +
  scale_fill_manual(name = "", values = c("90% interval" = "forestgreen")) +
  ylab("instanteneous growth rate (%)") +
  ylim(-5, 5)

gr_plot

combined_plot <- natural_plot / gr_plot + patchwork::plot_layout(guides = "collect",
  axes = "collect_x")

message("Saving outputs")


readr::write_csv(samples,
  file = here::here(
    "outputs", "data",
    "gt_samples.csv"
  )
)

readr::write_csv(gr_samples,
  file = here::here(
    "outputs", "data",
    "gt_gr_samples.csv"
  )
)

readr::write_csv(results,
  file = here::here(
    "outputs", "data",
    "gt_results.csv"
  )
)

readr::write_csv(gr_results,
  file = here::here(
    "outputs", "data",
    "gt_gr_results.csv"
  )
)

fs::dir_create(here::here(
  "outputs",
  "figures",
  "google_trends")
)

ggplot2::ggsave(
  filename = here::here(
    "outputs",
    "figures",
    "google_trends",
    "natural_gr_plot.png"
  ),
  plot = combined_plot,
  width = 8,
  height = 12)

message("GOOGLE TRENDS MODEL COMPLETE")
