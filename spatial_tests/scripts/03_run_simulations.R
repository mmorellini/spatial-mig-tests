# SCRIPT DESCRIPTION -----------------------------------------------------------
# 03_run_simulations.R
# Purpose:
#   - Use fitted PPML models (gravity and MST) for pred_q50
#   - Construct fitted values and 95% prediction intervals
#   - Run simulations of bilateral flows for each year
#   - Save gravity_year_sims and mst_year_sims for later use
#
# Inputs:
#   - data_processed/merged-flows.RData      (object: merged)
#   - data_processed/ppml_main_q50.RData     (objects: gravity_pois, mst_pois)
#   - scripts/utils/network_functions.R            (simulate_flows_years(), etc.)
#
# Outputs:
#   - data_processed/gravity_year_sims.rds
#   - data_processed/mst_year_sims.rds

library(dplyr)
library(fixest)
library(marginaleffects)
library(here)

# 0. Load helper functions, data, and models -----------------------------------

## Custom network / simulation functions
source(here::here("scripts", "utils", "network_functions.R"))

## Processed data (merged flows and covariates)
if (!exists("merged")) {
  load(here("data_processed", "merged-flows.RData"))
}

# Ensure a consistent dyad ordering within each year
merged <- merged %>%
  arrange(year, sending_country, receiving_country)

## PPML models estimated in Script 02_ppml_main_q50.R
## Objects expected in the .RData file: gravity_pois, mst_pois
load(here("data_processed", "ppml_main_q50.RData"))

# 1. Setup ---------------------------------------------------------------------

years <- sort(unique(merged$year))

country_names <- sort(unique(merged$sending_country))
n             <- length(country_names)

set.seed(415319416)   # For reproducibility
M_sim <- 10000         # Number of simulations per year 

conf_level_sim <- 0.95   # Must match conf_level used in predictions()

years <- sort(unique(merged$year))

# Build list of country pairs per year for simulations
# One data frame per year with dyads in the correct order
pairs_list <- lapply(years, function(y) {
  merged %>%
    filter(year == y) %>%
    select(sending_country, receiving_country)
})
names(pairs_list) <- as.character(years)

# 2. Gravity simulations -------------------------------------------------------

# 2.1 Fitted values and 95% prediction intervals
gravity_ci <- predictions(gravity_pois, conf_level = conf_level_sim)

merged <- merged %>%
  mutate(
    fitted_gravity   = gravity_ci$estimate,
    fitted_gravity_l = gravity_ci$conf.low,
    fitted_gravity_u = gravity_ci$conf.high
  )

# 2.2 Build per-year lists for simulate_flows_years()
fitted_gravity_list <- lapply(years, function(y) {
  merged$fitted_gravity[merged$year == y]
})
lower_gravity_list  <- lapply(years, function(y) {
  merged$fitted_gravity_l[merged$year == y]
})
upper_gravity_list  <- lapply(years, function(y) {
  merged$fitted_gravity_u[merged$year == y]
})
total_obs_list   <- lapply(years, function(y) {
  merged$pred_q50[merged$year == y]
})

names(fitted_gravity_list) <- as.character(years)
names(lower_gravity_list)  <- as.character(years)
names(upper_gravity_list)  <- as.character(years)
names(total_obs_list)   <- as.character(years)

# 2.3 Run simulations for Gravity
gravity_year_sims <- simulate_flows_years(
  fitted_values_list          = fitted_gravity_list,
  lower_bound_list            = lower_gravity_list,
  upper_bound_list            = upper_gravity_list,
  total_observed_flows_list   = total_obs_list,
  pairs_list                  = pairs_list,
  country_names               = country_names,
  years                       = years,
  M                           = M_sim
)

# Example access (optional check):
 gravity_year_sims[["2010"]][[1]]

# 3. MST simulations -----------------------------------------------------------

# 3.1 Fitted values and 95% prediction intervals
mst_ci <- predictions(mst_pois, conf_level = conf_level_sim)

merged <- merged %>%
  mutate(
    fitted_mst   = mst_ci$estimate,
    fitted_mst_l = mst_ci$conf.low,
    fitted_mst_u = mst_ci$conf.high
  )

# 3.2 Build per-year lists for MST
fitted_mst_list <- lapply(years, function(y) {
  merged$fitted_mst[merged$year == y]
})
lower_mst_list  <- lapply(years, function(y) {
  merged$fitted_mst_l[merged$year == y]
})
upper_mst_list  <- lapply(years, function(y) {
  merged$fitted_mst_u[merged$year == y]
})

names(fitted_mst_list) <- as.character(years)
names(lower_mst_list)  <- as.character(years)
names(upper_mst_list)  <- as.character(years)

# 3.3 Run simulations for MST
mst_year_sims <- simulate_flows_years(
  fitted_values_list          = fitted_mst_list,
  lower_bound_list            = lower_mst_list,
  upper_bound_list            = upper_mst_list,
  total_observed_flows_list   = total_obs_list,
  pairs_list                  = pairs_list,
  country_names               = country_names,
  years                       = years,
  M                           = M_sim
)

# Example access (optional check):
#mst_year_sims[["2010"]][[1]]

# 4. Save simulation objects ---------------------------------------------------

if (!dir.exists(here("data_processed"))) {
  dir.create(here("data_processed"), recursive = TRUE)
}

saveRDS(gravity_year_sims, here("data_processed", "gravity_year_sims.rds"))
saveRDS(mst_year_sims,     here("data_processed", "mst_year_sims.rds"))
