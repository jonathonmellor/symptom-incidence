# Symptom Incidence Google WCIS

Code for the upcoming manuscript: "Evaluating Google Trends as a Proxy for Symptom Incidence: Insights from the Winter COVID Infection Study in England 2023/24"

Authors: Phoebe Asplin, Martyn Fyles Jack Kennedy, Thomas Ward, Jonathon Mellor.

Due to the disclosive nature of the patient level data used in this project the data is only available upon request. The code in this repository therefore does not run and is instead provided by the authors to illustrate the processing and modelling completed.

# Structure

## WCIS

Data are cleaned in `./load_data.R`.

To calculate incidence we need both prevalence and symptom duration estimates.

`./prevalence.R` and `./duration.R` both run and store csv files with model samples for each symptom and subgroup.

`./incidence.R` combines the prevalence and duration estimates to give an incidence by symptom.

## Google Symptoms

Data are loaded and modelled within `./google_trends.R`, applied to each symptom with model samples stored to csv.

## Analysis

Each model run script has some diagnostic output.

The primary analysis happens in `./compare.R` (which runs all upstream models) where the samples from the WCIS incidence model are compared against the Google GAM.
This includes plotting both data against each other on natural and growth rate scales, as well as cross correlation analysis.

Descriptive analysis is contained within `./cohort_summary.R` to demonstrate the WCIS cohort.

## Dependency

There are multiple models feeding into the final analysis.

The WCIS incidence model (`incidence.R`) depends on `duration.R` and `prevalence.R`.

The Google incidence model relies only on `google_trends.R`.

The comparison of both (`compare.R`) relies on `incidence.R` and `google_trends.R`.

