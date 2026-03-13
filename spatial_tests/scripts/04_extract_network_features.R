# SCRIPT DESCRIPTION -----------------------------------------------------------
# 04_extract_network_features.R
# Purpose:
#   - Build observed flow matrices by year from `merged` (pred_q50).
#   - Compute global network metrics (reciprocity, Gini, ACV, migration inequality)
#     for each year, simulation, and model (Gravity / MST), using parallelism.
#   - Additionally compute the same metrics for observed flows at multiple
#     Bayesian quantiles (eg pred_q25, pred_q50, pred_q75).
#
# Inputs:
#   - data_processed/merged-flows.RData        (object: merged)
#   - data_processed/gravity_year_sims.rds     (list: year -> list of sim matrices)
#   - data_processed/mst_year_sims.rds         (list: year -> list of sim matrices)
#   - scripts/utils/network_functions.R        (country_names, create_flow_matrix,
#                                               reciprocity_min)
# Outputs:
#   - data_processed/global_values_long.rds
#   - data_processed/observed_matrices_year.rds        (pred_q50 matrices)
#   - data_processed/observed_values_all_q.rds  (indices for all pred_q*)

library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)  # for parallel processing
library(rlist)
library(migration.indices) # for migration.gini.total, migration.acv, migration.inequality

# 0. Helper functions ----------------------------------------------------------

source(here::here("scripts", "utils", "network_functions.R"))
# Expecting: country_names (optional), create_flow_matrix(), reciprocity_min()

# 1. Load data and simulations -------------------------------------------------

if (!exists("merged")) {
  load(here("data_processed", "merged-flows.RData"))
}

gravity_year_sims <- readRDS(here("data_processed", "gravity_year_sims.rds"))
mst_year_sims     <- readRDS(here("data_processed", "mst_year_sims.rds"))

# Ensure common set of years
years_chr <- intersect(names(gravity_year_sims), names(mst_year_sims))
years_chr <- sort(years_chr)
years     <- as.integer(years_chr)

# Define number and names of countries
country_names <- sort(unique(merged$sending_country))
n             <- length(country_names)

# 2. Build observed matrices by year from `merged` (pred_q50 only) ------------

obs_matrices_year <- map(years, function(y) {
  flows_y_df <- merged %>%
    filter(year == y) %>%
    select(sending_country, receiving_country, pred_q50)
  
  create_flow_matrix(
    flows_df      = flows_y_df,
    country_names = country_names,
    value_col     = "pred_q50"
  )
})

names(obs_matrices_year) <- years_chr

# 3. Global metrics via parallel processing (simulations) ----------------------

plan(multisession)

gravity_metrics_list <- future_map(gravity_year_sims, compute_metrics, .progress = TRUE)
mst_metrics_list     <- future_map(mst_year_sims,     compute_metrics, .progress = TRUE)

plan(sequential)

gravity_long <- combine_metrics(gravity_metrics_list, "Gravity")
mst_long     <- combine_metrics(mst_metrics_list,     "MST")

global_values_long <- bind_rows(gravity_long, mst_long)

global_values_long <- global_values_long %>%
  mutate(Year = as.character(Year))

# 4. Observed global values by year, for multiple pred_q* ----------------------
#    - Look for columns named pred_q25, pred_q50, pred_q75 (or any pred_q*) and
#      compute the same four indices for each year × quantile.

# Which pred_q* columns are available?
pred_cols_all <- grep("^pred_q", names(merged), value = TRUE)

# Keeping only the main three:
# pred_cols_all <- intersect(c("pred_q25", "pred_q50", "pred_q75"), names(merged))

if (length(pred_cols_all) == 0L) {
  stop("No pred_q* columns found in `merged`. Cannot compute observed indices by quantile.")
}

# Compute indices for each quantile & year
obs_values_all_q <- map_dfr(pred_cols_all, function(vcol) {
  map_dfr(years_chr, function(yr_chr) {
    
    flows_y_df <- merged %>%
      filter(year == as.integer(yr_chr)) %>%
      select(sending_country, receiving_country, all_of(vcol))
    
    mat <- create_flow_matrix(
      flows_df      = flows_y_df,
      country_names = country_names,
      value_col     = vcol
    )
    
    tibble(
      Year         = yr_chr,
      quantile_var = vcol,
      # nicer label, e.g. "q25", "q50", "q75"
      quantile     = sub("^pred_", "", vcol),
      statistic_raw = c("gini", "acv", "inequality", "reciprocity"),
      obs_value     = c(
        migration.gini.total(mat),
        migration.acv(mat),
        migration.inequality(mat) / sum(rowSums(mat)),
        reciprocity_min(mat)
      )
    )
  })
})

# For compatibility: keep only pred_q50 when joining to sims
obs_values_by_year <- obs_values_all_q %>%
  filter(quantile_var == "pred_q50") %>%
  select(Year, statistic_raw, obs_value)

# 5. Add means per model & year, and merge with observed (q50) -----------------

dummy <- global_values_long %>%
  group_by(statistic_raw, Model, Year) %>%
  summarize(mean = mean(value), .groups = "drop") %>%
  pivot_wider(
    names_from  = Model,
    values_from = mean,
    names_prefix = "mean_"
  ) %>%
  left_join(obs_values_by_year, by = c("Year", "statistic_raw"))

global_values_long <- global_values_long %>%
  left_join(dummy, by = c("Year", "statistic_raw"))

# Rename statistic nicely (Gini, ACV, Migration inequality, Reciprocity)
global_values_long <- global_values_long %>%
  mutate(
    statistic = recode(
      statistic_raw,
      gini        = "Gini",
      acv         = "ACV",
      inequality  = "Migration inequality",
      reciprocity = "Reciprocity"
    ),
    statistic = factor(
      statistic,
      levels = c("Gini", "ACV", "Migration inequality", "Reciprocity")
    )
  )

# 6. Save outputs --------------------------------------------------------------

if (!dir.exists(here("data_processed"))) {
  dir.create(here("data_processed"), recursive = TRUE)
}

saveRDS(global_values_long,
        here("data_processed", "global_values_long.rds"))

# Pred_q50 matrices (used by other scripts)
saveRDS(obs_matrices_year,
        here("data_processed", "observed_matrices_year.rds"))

# New: all quantiles’ observed indices
saveRDS(obs_values_all_q,
        here("data_processed", "observed_values_all_q.rds"))