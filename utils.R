# storage of utility artifacts


#  age classification ####

age_groups_mutate <- function(df) { # Use this function to define the age groups we want to model
  df <- df |> mutate(age_group = case_when(
    # ONS proposed age bands
    age < 18 ~ "3 to 17 years",
    age <= 34 ~ "18 to 34 years",
    age <= 44 ~ "35 to 44 years",
    age <= 54 ~ "45 to 54 years",
    age <= 64 ~ "55 to 64 years",
    age <= 74 ~ "65 to 74 years",
    TRUE ~ "75 years and over"
  ))
  return(df)
}

# Generated ordered vector of region names
# probably better as a config
region_order <- c(
  "North East",
  "North West",
  "Yorkshire and The Humber",
  "East Midlands",
  "West Midlands",
  "East of England",
  "London",
  "South East",
  "South West"
)

# used to arrange ages
age_group_order <- c(
  "3 to 17 years",
  "18 to 34 years",
  "35 to 44 years",
  "45 to 54 years",
  "55 to 64 years",
  "65 to 74 years",
  "75 years and over"
)
