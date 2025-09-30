# compare WCIS symptom incidence vs google trends data

# depends on:
# - incidence.R
# - google_trends.R

pacman::p_load(dplyr, ggh4x)

# note: takes ~15 mins
source(here::here("incidence.R"))
# note: takes ~30s
source(here::here("google_trends.R"))


source(here::here("utils.R"))
data_order <- c("Google trends", "All", age_group_order)

ggplot2::theme_set(ggplot2::theme_bw())

fs::dir_create(here::here(
  "outputs",
  "figures",
  "compare"
))

fs::dir_create(here::here(
  "outputs",
  "tables"
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


## LOAD DATA ##

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
  dplyr::mutate(data = age_group)


# WCIS growth rates
wcis_gr_samples <- wcis_samples |>
  dplyr::arrange(symptom, .draw, age_group, date) |>
  dplyr::mutate(
    diff_time = as.numeric(date - lag(date)),
    diff_incidence = incidence - lag(incidence),
    wcis_gr = (diff_incidence / diff_time) / incidence * 100
  ) |>
  dplyr::group_by(symptom, age_group, .draw) |>
  dplyr::filter(date > min(date)) |>
  dplyr::ungroup()

# calculate quantiles
wcis_gr_summary <- wcis_gr_samples |>
  dplyr::summarise(
    pi_50 = quantile(wcis_gr, 0.5),
    pi_95 = quantile(wcis_gr, 0.95),
    pi_5 = quantile(wcis_gr, 0.05),
    pi_75 = quantile(wcis_gr, 0.75),
    pi_25 = quantile(wcis_gr, 0.25),
    .by = c("symptom", "age_group", "date")
  ) |>
  dplyr::mutate(data = age_group)


# load google trends
gt_samples <- readr::read_csv(
  file = here::here(
    "outputs",
    "data",
    "gt_samples.csv"
  )
) |>
  dplyr::mutate(
    data = "Google trends",
    symptom = replace(symptom, symptom == "shortness_of_breath", "short_breath"),
    symptom = replace(symptom, symptom == "anosmia", "loss_smell")
  ) |>
  # filter to only selected symptoms
  dplyr::filter(stringr::str_detect(symptom, paste(symptom_opts, collapse = "|")))

# calculate quantiles
gt_summary <- gt_samples |>
  dplyr::summarise(
    pi_50 = quantile(.fitted, 0.5),
    pi_95 = quantile(.fitted, 0.95),
    pi_5 = quantile(.fitted, 0.05),
    pi_75 = quantile(.fitted, 0.75),
    pi_25 = quantile(.fitted, 0.25),
    value = unique(value),
    .by = c("date", "symptom")
  ) |>
  dplyr::mutate(data = "Google trends")


# load google trends growth rates
gt_gr_samples <- gt_samples |>
  dplyr::arrange(symptom, .draw, date) |>
  dplyr::mutate(
    diff_time = as.numeric(date - lag(date)),
    diff_incidence = .fitted - lag(.fitted),
    gt_gr = (diff_incidence / diff_time) / .fitted * 100
  ) |>
  dplyr::group_by(symptom, .draw) |>
  dplyr::filter(date > min(date)) |>
  dplyr::ungroup()

# calculate quantiles
gt_gr_summary <- gt_gr_samples |>
  dplyr::summarise(
    pi_50 = quantile(gt_gr, 0.5),
    pi_95 = quantile(gt_gr, 0.95),
    pi_5 = quantile(gt_gr, 0.05),
    pi_75 = quantile(gt_gr, 0.75),
    pi_25 = quantile(gt_gr, 0.25),
    .by = c("date", "symptom")
  ) |>
  dplyr::mutate(data = "Google trends")


## PLOTS ##

# plot incidence vs google trends
incidence_plot <- dplyr::bind_rows(wcis_summary, gt_summary) |>
  ggplot() +
  geom_ribbon(aes(x = date, ymax = pi_95, ymin = pi_5, fill = symptom, group = symptom),
    alpha = 0.6
  ) +
  geom_line(aes(x = date, y = pi_50, group = symptom)) +
  labs(x = "Date", y = "Incidence rate of symptom per day") +
  theme(axis.text.x = element_text(size = 8, angle = 90)) +
  scale_x_date(
    breaks = "3 weeks", minor_breaks = "1 week",
    labels = scales::date_format("%b-%d")
  ) +
  scale_y_continuous(labels = scales::comma) +
  ggh4x::facet_grid2(
    factor(data, data_order) ~ symptom,
    scales = "free_y",
    independent = "y",
    labeller = labeller(symptom = as_labeller(symptom_labels))
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

incidence_plot

ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "compare",
    "incidence_comparison.png"
  ),
  incidence_plot,
  height = 14,
  width = 12
)


# plot growth rates
gr_plot <- dplyr::bind_rows(wcis_gr_summary, gt_gr_summary) |>
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_ribbon(aes(x = date, ymax = pi_95, ymin = pi_5, fill = symptom, group = symptom),
    alpha = 0.6
  ) +
  geom_line(aes(x = date, y = pi_50, group = symptom)) +
  labs(x = "Date", y = "Growth rate of symptom incidence per day") +
  theme(axis.text.x = element_text(size = 8, angle = 90)) +
  scale_x_date(
    breaks = "3 weeks", minor_breaks = "1 week",
    labels = scales::date_format("%b-%d")
  ) +
  scale_y_continuous(labels = scales::comma) +
  ggh4x::facet_grid2(
    factor(data, data_order) ~ symptom,
    scales = "free_y",
    independent = "y",
    labeller = labeller(symptom = as_labeller(symptom_labels))
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

gr_plot

ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "compare",
    "gr_comparison.png"
  ),
  gr_plot,
  height = 14,
  width = 12
)


# ccf for incidence
ccf_incidence <- dplyr::inner_join(wcis_samples, gt_samples, by = c("symptom", "date", ".draw")) |>
  dplyr::arrange(symptom, age_group, date) |>
  dplyr::group_by(symptom, age_group, .draw) |>
  dplyr::reframe(
    ccf = c(ccf(incidence, .fitted, plot = FALSE))$acf,
    lag = c(ccf(incidence, .fitted, plot = FALSE))$lag
  ) |>
  dplyr::ungroup()

ccf_incidence_plot <- ccf_incidence |>
  dplyr::summarise(
    pi_50 = quantile(ccf, 0.5),
    pi_95 = quantile(ccf, 0.95),
    pi_5 = quantile(ccf, 0.05),
    pi_75 = quantile(ccf, 0.75),
    pi_25 = quantile(ccf, 0.25),
    .by = c("symptom", "age_group", "lag")
  ) |>
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_ribbon(aes(x = lag, ymin = pi_5, ymax = pi_95, fill = symptom), alpha = 0.6) +
  geom_line(aes(x = lag, y = pi_50)) +
  labs(y = "Cross correlation of WCIS incidence rate and Google trends", x = "Lag (days)") +
  scale_x_continuous(breaks = seq(-15, 15, 5), minor_breaks = NULL) +
  ylim(-1, 1) +
  facet_grid(
    factor(age_group, c("All", age_group_order)) ~ symptom,
    labeller = labeller(symptom = as_labeller(symptom_labels))
  ) +
  theme(legend.position = "none")

ccf_incidence_plot

ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "compare",
    "incidence_ccf.png"
  ),
  ccf_incidence_plot,
  height = 12,
  width = 12
)

sf_digits <- 2
sf_format <- "fg"

max_lags_incidence <- ccf_incidence |>
  dplyr::summarise(
    lag = lag[which.max(ccf)],
    ccf = max(ccf),
    .by = c("symptom", "age_group", ".draw")
  ) |>
  dplyr:::summarise(
    lag_pi_50 = quantile(lag, 0.5),
    lag_pi_95 = quantile(lag, 0.95),
    lag_pi_5 = quantile(lag, 0.05),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::summarise(
    value =
      glue::glue(
        "{formatC(signif(lag_pi_50, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')}",
        " [{formatC(signif(lag_pi_5, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')},",
        " {formatC(signif(lag_pi_95, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')}]"
      ),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::mutate(age_group = factor(age_group, levels = c("All", age_group_order))) |>
  dplyr::arrange(age_group) |>
  tidyr::pivot_wider(
    names_from = age_group,
    values_from = value
  )

max_lags_incidence


readr::write_csv(max_lags_incidence,
  file = here::here(
    "outputs",
    "tables",
    "max_lags_incidence.csv"
  )
)

max_corrs_incidence <- ccf_incidence |>
  dplyr::summarise(
    lag = lag[which.max(ccf)],
    ccf = max(ccf),
    .by = c("symptom", "age_group", ".draw")
  ) |>
  dplyr:::summarise(
    ccf_pi_50 = quantile(ccf, 0.5),
    ccf_pi_95 = quantile(ccf, 0.95),
    ccf_pi_5 = quantile(ccf, 0.05),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::summarise(
    value = glue::glue(
      "{formatC(signif(ccf_pi_50, digits=sf_digits)",
      ", format=sf_format, digits=sf_digits,flag='#')}",
      " [{formatC(signif(ccf_pi_5, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')},",
      " {formatC(signif(ccf_pi_95, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')}]"
    ),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::mutate(age_group = factor(age_group, levels = c("All", age_group_order))) |>
  dplyr::arrange(age_group) |>
  tidyr::pivot_wider(
    names_from = age_group,
    values_from = value
  )

max_corrs_incidence

readr::write_csv(max_corrs_incidence,
  file = here::here(
    "outputs",
    "tables",
    "max_corrs_incidence.csv"
  )
)



# ccf for growth rate
ccf_gr <- dplyr::inner_join(wcis_gr_samples, gt_gr_samples, by = c("symptom", "date", ".draw")) |>
  dplyr::arrange(symptom, age_group, date) |>
  dplyr::group_by(symptom, age_group, .draw) |>
  dplyr::reframe(
    ccf = c(ccf(wcis_gr, gt_gr, plot = FALSE))$acf,
    lag = c(ccf(wcis_gr, gt_gr, plot = FALSE))$lag
  ) |>
  dplyr::ungroup()

ccf_gr_plot <- ccf_gr |>
  dplyr::summarise(
    pi_50 = quantile(ccf, 0.5),
    pi_95 = quantile(ccf, 0.95),
    pi_5 = quantile(ccf, 0.05),
    pi_75 = quantile(ccf, 0.75),
    pi_25 = quantile(ccf, 0.25),
    .by = c("symptom", "age_group", "lag")
  ) |>
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_ribbon(aes(x = lag, ymin = pi_5, ymax = pi_95, fill = symptom), alpha = 0.6) +
  geom_line(aes(x = lag, y = pi_50)) +
  labs(y = "Cross correlation of WCIS growth rate vs Google trends", x = "Lag (days)") +
  scale_x_continuous(breaks = seq(-15, 15, 5), minor_breaks = NULL) +
  ylim(-1, 1) +
  facet_grid(
    factor(age_group, c("All", age_group_order)) ~ symptom,
    labeller = labeller(symptom = as_labeller(symptom_labels))
  ) +
  theme(legend.position = "none")

ccf_gr_plot

ggplot2::ggsave(
  here::here(
    "outputs",
    "figures",
    "compare",
    "gr_ccf.png"
  ),
  ccf_gr_plot,
  height = 12,
  width = 12
)

max_lags_gr <- ccf_gr |>
  dplyr::summarise(
    lag = lag[which.max(ccf)],
    ccf = max(ccf),
    .by = c("symptom", "age_group", ".draw")
  ) |>
  dplyr:::summarise(
    lag_pi_50 = quantile(lag, 0.5),
    lag_pi_95 = quantile(lag, 0.95),
    lag_pi_5 = quantile(lag, 0.05),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::summarise(
    value = glue::glue(
      "{formatC(signif(lag_pi_50, digits=sf_digits),",
      " format=sf_format, digits=sf_digits,flag='#')}",
      " [{formatC(signif(lag_pi_5, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')},",
      " {formatC(signif(lag_pi_95, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')}]"
    ),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::mutate(age_group = factor(age_group, levels = c("All", age_group_order))) |>
  dplyr::arrange(age_group) |>
  tidyr::pivot_wider(
    names_from = age_group,
    values_from = value
  )

max_lags_gr

readr::write_csv(max_lags_gr,
  file = here::here(
    "outputs",
    "tables",
    "max_lags_gr.csv"
  )
)

max_corrs_gr <- ccf_gr |>
  dplyr::summarise(
    lag = lag[which.max(ccf)],
    ccf = max(ccf),
    .by = c("symptom", "age_group", ".draw")
  ) |>
  dplyr:::summarise(
    ccf_pi_50 = quantile(ccf, 0.5),
    ccf_pi_95 = quantile(ccf, 0.95),
    ccf_pi_5 = quantile(ccf, 0.05),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::summarise(
    value = glue::glue(
      "{formatC(signif(ccf_pi_50, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')}",
      " [{formatC(signif(ccf_pi_5, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')},",
      " {formatC(signif(ccf_pi_95, digits=sf_digits), format=sf_format, digits=sf_digits,flag='#')}]"
    ),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::mutate(age_group = factor(age_group, levels = c("All", age_group_order))) |>
  dplyr::arrange(age_group) |>
  tidyr::pivot_wider(
    names_from = age_group,
    values_from = value
  )

readr::write_csv(max_corrs_gr,
  file = here::here(
    "outputs",
    "tables",
    "max_corrs_gr.csv"
  )
)

## Explore max correlations and max lags in 2D ----
# Not used in the paper.

combined_gr_ccf <- ccf_gr |>
  dplyr::summarise(
    lag = lag[which.max(ccf)],
    ccf = max(ccf),
    .by = c("symptom", "age_group", ".draw")
  )

combined_gr_ccf_summarised <- combined_gr_ccf |>
  dplyr:::summarise(
    ccf_pi_50 = quantile(ccf, 0.5),
    ccf_pi_95 = quantile(ccf, 0.95),
    ccf_pi_5 = quantile(ccf, 0.05),
    lag_pi_50 = quantile(lag, 0.5),
    lag_pi_95 = quantile(lag, 0.95),
    lag_pi_5 = quantile(lag, 0.05),
    .by = c("symptom", "age_group")
  ) |>
  dplyr::mutate(
    age_group = factor(age_group, c("All", age_group_order)),
    symptom = factor(symptom, labels = symptom_labels)
  )

colours <- c("black", RColorBrewer::brewer.pal(n = 9, "Spectral"))
colours <- colours[-c(5, 6)]
colours[5] <- "#BA8E23"

jitter <- position_jitter(
  seed = 123,
  width = 0.5,
  height = 0.05
)

# facet by symptom to allow age group comparison within a symptom group
combined_gr_ccf_summarised |>
  dplyr::mutate(is_national = dplyr::if_else(age_group == "All", "A", "B")) |>
  ggplot(aes(x = lag_pi_50, y = ccf_pi_50)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha = 0.4) +
  geom_linerange(aes(ymin = ccf_pi_5, ymax = ccf_pi_95, color = age_group), alpha = 1, position = jitter) +
  geom_linerange(aes(xmin = lag_pi_5, xmax = lag_pi_95, color = age_group), alpha = 1, position = jitter) +
  geom_point(aes(color = age_group, size = is_national), position = jitter) +
  scale_x_continuous(breaks = seq(-15, 15, 5), minor_breaks = seq(-20, 20, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_size_manual(values = c("A" = 2, "B" = 1.5)) +
  scale_shape_manual(values = c("A" = 19, "B" = 1)) +
  labs(
    y = "Maximum cross correlation",
    x = "Lag (days) of maximum cross correlation"
  ) +
  scale_color_manual(name = "Age group", values = colours) +
  facet_wrap(~symptom) +
  theme(legend.position = "bottom") +
  guides(
    size = "none",
    shape = "none",
    color = guide_legend(nrow = 2, byrow = TRUE)
  )

message("COMPARISON COMPLETE")
