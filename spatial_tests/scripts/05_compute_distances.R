# SCRIPT DESCRIPTION -----------------------------------------------------------
# 05_compute_distances.R
# Purpose:
#   - Compute distance-based and other diagnostics between observed and simulated flows:
#       1) By country over time (inflows) for each model (Mahalanobis distance)
#       2) For all corridors over time for each model (various fit statistics)
#   - Save tidy outputs for ggplot-based visualisation.
#
# Inputs:
#   - data_processed/gravity_year_sims.rds
#   - data_processed/mst_year_sims.rds
#   - data_processed/observed_matrices_year.rds
#   - data_processed/merged-flows.RData
#   - scripts/utils/network_functions.R
#   - scripts/utils/compute_mahalanobis.R
#
# Outputs:
#   - data_processed/country_distances_long.rds
#   - data_processed/sim_summaries.rds

library(here)
library(dplyr)
library(purrr)
library(tidyr)
library(countrycode)

# 0. Source helper scripts -----------------------------------------------------

source(here("scripts", "utils", "network_functions.R"))
source(here("scripts", "utils", "compute_mahalanobis.R"))

# 1. Load data and simulations -------------------------------------------------

# Observed matrices by year
if (!exists("obs_matrices_year")) {
  obs_matrices_year <- readRDS(
    here("data_processed", "observed_matrices_year.rds")
  )
}

# Migration flows (Bayesian)
if (!exists("merged")) {
  load(here("data_processed", "merged-flows.RData"))  # loads `merged`
}

# Simulations by year
gravity_year_sims <- readRDS(here("data_processed", "gravity_year_sims.rds"))
mst_year_sims     <- readRDS(here("data_processed", "mst_year_sims.rds"))

# Ensure common years (character labels)
years_chr <- intersect(names(gravity_year_sims), names(mst_year_sims))
years_chr <- sort(years_chr)
years     <- as.integer(years_chr)

# 2. By country: Mahalanobis over time -----------------------------------------
#    Using get_mahalanobis_by_country():
#      - obs: list(year -> observed matrix)
#      - sim: list(year -> list of simulated matrices)
#    Returns:
#      year, Year, country, model, mahalanobis distance

mhd_country_gravity <- get_mahalanobis_by_country(
  obs        = obs_matrices_year,
  sim        = gravity_year_sims,
  model_name = "Gravity",
  verbose    = TRUE
)   

mhd_country_mst <- get_mahalanobis_by_country(
  obs        = obs_matrices_year,
  sim        = mst_year_sims,
  model_name = "MST",
  verbose    = TRUE
)

mhd_country_gravity <- mhd_country_gravity %>%
  dplyr::group_by(year, country, model) %>%
  dplyr::summarise(
    mahalanobis = mean(mahalanobis, na.rm = TRUE),
    .groups = "drop"
  )

mhd_country_mst <- mhd_country_mst %>%
  dplyr::group_by(year, country, model) %>%
  dplyr::summarise(
    mahalanobis = mean(mahalanobis, na.rm = TRUE),
    .groups = "drop"
  )

country_distances_long <- bind_rows(
  mhd_country_gravity,
  mhd_country_mst) %>%
  mutate(
    Year = as.numeric(year)) %>%
  relocate(Year, .before = year)

# 3. By corridor: Fit statistics over time -------------------------------------

# Get summary statistics by corridor across simulations
summaries_gravity <- map_dfr(names(gravity_year_sims), function(yr_chr) {
  sims_yr <- gravity_year_sims[[yr_chr]]  # list of matrices
  mat_first <- sims_yr[[1]]
  origin_codes <- rownames(mat_first)
  dest_codes   <- colnames(mat_first)
  
  expand_grid(
    origin = origin_codes,
    destination = dest_codes
  ) %>%
    # drop diagonals if you want
    filter(origin != destination) %>%
    rowwise() %>%
    mutate(
      # extract all sims for this (origin, destination, year)
      values = list(vapply(
        sims_yr,
        function(M) M[origin, destination],
        numeric(1)
      )),
      mean_sim = mean(values),
      sd_sim   = sd(values),
      min_sim  = min(values),
      max_sim  = max(values),
      q05      = quantile(values, 0.05),
      q25      = quantile(values, 0.25),
      q50      = quantile(values, 0.50),
      q75      = quantile(values, 0.75),
      q95      = quantile(values, 0.95),
      min_sim  = min(values),
      max_sim  = max(values)
    ) %>%
    ungroup() %>%
    select(-values) %>%
    mutate(
      model = "Gravity",
      year  = as.integer(yr_chr)
    )
})

summaries_mst <- map_dfr(names(mst_year_sims), function(yr_chr) {
  sims_yr <- mst_year_sims[[yr_chr]]  # list of matrices
  mat_first <- sims_yr[[1]]
  origin_codes <- rownames(mat_first)
  dest_codes   <- colnames(mat_first)
  
  expand_grid(
    origin = origin_codes,
    destination = dest_codes
  ) %>%
    # drop diagonals
    filter(origin != destination) %>%
    rowwise() %>%
    mutate(
      values = list(vapply(
        sims_yr,
        function(M) M[origin, destination],
        numeric(1)
      )),
      mean_sim = mean(values),
      sd_sim   = sd(values),
      min_sim  = min(values),
      max_sim  = max(values),
      q05      = quantile(values, 0.05),
      q25      = quantile(values, 0.25),
      q50      = quantile(values, 0.50),
      q75      = quantile(values, 0.75),
      q95      = quantile(values, 0.95),
      min_sim  = min(values),
      max_sim  = max(values)
    ) %>%
    ungroup() %>%
    select(-values) %>%
    mutate(
      model = "MST",
      year  = as.integer(yr_chr)
    )
})

sim_summaries <- bind_rows(summaries_gravity, summaries_mst)

# Add regional information

selected_corridors <- countrycode::codelist %>%
  filter(!is.na(iso2c) & !is.na(un.regionsub.name)) %>%
  select(iso3c, iso2c, un.regionsub.name)

selected_corridors$iso2c[selected_corridors$iso2c == "GB"] <- "UK"  # UK countrycode fix

selected_corridors$un.regionsub.name[selected_corridors$iso2c == "CY"] <- "Southern Europe"  # region name fix


corridors_reg <- merged %>%
  select(iso3c_o, iso3c_d) %>% 
  distinct() %>% 
  # Origin region
  left_join(selected_corridors,
            by = c("iso3c_o" = "iso3c")) %>%
  rename(origin_region = un.regionsub.name,
         origin = iso2c) %>% 
  # Destination region
  left_join(selected_corridors,
            by = c("iso3c_d" = "iso3c")) %>%
  rename(dest_region = un.regionsub.name,
         destination = iso2c) %>% 
  # region-to-region corridor label
  mutate(
    corridor_pair = paste(
      origin,
      destination,
      sep = " \u2192 "
    ),
    region_pair = paste(
      origin_region,
      dest_region,
      sep = " \u2192 "
    )
  )


sim_summaries <- sim_summaries %>%
  left_join(corridors_reg) %>% 
  left_join(merged %>% 
              select(year, iso3c_o, iso3c_d, starts_with("pred_")))

corridor_fit <- sim_summaries %>%
  mutate(
    abs_err = abs(mean_sim - pred_q50),
    z_score = if_else(
      !is.na(sd_sim) & sd_sim > 0,
      (pred_q50 - mean_sim) / sd_sim,
      NA_real_
    ),
    # Approximate rank probability: P(sim ≤ obs) under Normal(mean_sim, sd_sim)
    rank_prob = if_else(
      !is.na(sd_sim) & sd_sim > 0,
      pnorm(pred_q50, mean = mean_sim, sd = sd_sim),
      NA_real_
    )
  )

# Add regional abbreviations:

region_abbrev <- c(
  "Eastern Europe"  = "EE",
  "Northern Europe" = "NE",
  "Southern Europe" = "SE",
  "Western Europe"  = "WE"
)

corridor_fit <- corridor_fit %>%
  mutate(
    region_short_o = region_abbrev[origin_region],
    region_short_d = region_abbrev[dest_region],
    region_corr    = paste(region_short_o, "→", region_short_d)
  )

corridor_features <- corridor_fit %>%
  group_by(year, origin, destination, model) %>%
  summarise(
    mean_abs_err   = mean(abs_err, na.rm = TRUE),
    mean_abs_z     = mean(abs(z_score), na.rm = TRUE),
    mean_rank_prob = mean(rank_prob,  na.rm = TRUE),
    .groups = "drop"
  )

corridor_features <- corridor_features %>%
  pivot_wider(
    names_from  = model,
    values_from = c(mean_abs_err, mean_abs_z, mean_rank_prob
    ),
    names_sep   = "_"
  ) %>% 
  left_join(sim_summaries %>% 
              select(year, origin, destination, corridor_pair, region_pair))

corridor_model_comp <- corridor_features %>%
  mutate(
    # smaller error is better
    better_model = case_when(
      is.na(mean_abs_err_Gravity) | is.na(mean_abs_err_MST) ~ NA_character_,
      mean_abs_err_Gravity + mean_abs_err_MST == 0          ~ NA_character_,
      mean_abs_err_Gravity < mean_abs_err_MST               ~ "Gravity",
      mean_abs_err_MST    < mean_abs_err_Gravity            ~ "MST",
      TRUE                                                   ~ "Tie"
    ),
    diff_abs_err = mean_abs_err_MST - mean_abs_err_Gravity  # >0 => Gravity better
  )

# Sanity check: model performance over corridors
check <- corridor_model_comp %>%
  count(origin, destination, better_model) %>%
  group_by(origin, destination) %>%
  mutate(pct = n / sum(n)) %>%
  arrange(origin, destination, desc(pct))

# 5. Save outputs --------------------------------------------------------------

if (!dir.exists(here("data_processed"))) {
  dir.create(here("data_processed"), recursive = TRUE)
}

saveRDS(
  country_distances_long,
  here("data_processed", "country_distances_long.rds")
)

saveRDS(
  sim_summaries,
  here("data_processed", "sim_summaries.rds")
)

saveRDS(
  corridor_fit,
  here("data_processed", "corridor_fit.rds")
)

