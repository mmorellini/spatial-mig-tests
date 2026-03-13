# network_functions.R
# Purpose: helper functions for simulations and network metrics.

library(truncnorm)

# -------------------------------------------------------------------
# 1. Countries and matrix constructors
# -------------------------------------------------------------------

country_names <- c(
  "AT", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES",
  "FI", "FR", "GR", "HR", "HU", "IE", "IS", "IT", "LT", "LU",
  "LV", "MT", "NL", "NO", "PL", "PT", "RO", "SE", "SI", "SK", "UK"
)

n <- length(country_names)

# Generic function to create a mobility / flow matrix from a vector
create_flow_matrix <- function(flows_df, country_names, value_col) {
  n <- length(country_names)
  
  # Zero matrix with row/column names
  mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(country_names, country_names)
  )
  
  # Map dyads to matrix holes
  idx_send <- match(flows_df$sending_country, country_names)
  idx_recv <- match(flows_df$receiving_country, country_names)
  
  mat[cbind(idx_send, idx_recv)] <- flows_df[[value_col]]
  mat
}

# -------------------------------------------------------------------
# 2. Standard deviation from confidence interval
# -------------------------------------------------------------------

calculate_sd <- function(lower_bound, upper_bound, confidence_level = 0.95) {
  # Generic z from normal
  z_value <- qnorm((1 + confidence_level) / 2)
  
  if (any(is.na(z_value))) {
    stop("Unsupported confidence level: ", confidence_level)
  }
  
  (upper_bound - lower_bound) / (2 * z_value)
}

# -------------------------------------------------------------------
# 3. Simulate flows over years
# -------------------------------------------------------------------
simulate_flows_years <- function(fitted_values_list,
                                 lower_bound_list,
                                 upper_bound_list,
                                 total_observed_flows_list,
                                 pairs_list,
                                 country_names,
                                 years,
                                 M = 10000) {
  all_year_sims <- vector("list", length(years))
  names(all_year_sims) <- as.character(years)
  
  for (i in seq_along(years)) {
    yr      <- years[i]
    yr_chr  <- as.character(yr)
    
    fitted_values        <- fitted_values_list[[yr_chr]]
    lower_bound          <- lower_bound_list[[yr_chr]]
    upper_bound          <- upper_bound_list[[yr_chr]]
    total_observed_flows <- total_observed_flows_list[[yr_chr]]
    pairs_df             <- pairs_list[[yr_chr]]   # must match the same order
    
    # Pre-compute SD vector (same for all simulations in year yr)
    sd_vec <- calculate_sd(
      lower_bound,
      upper_bound,
      confidence_level = 0.95
    )
    
    year_sims <- vector("list", M)
    
    for (m in seq_len(M)) {
      # 1) Simulate flows (vector)
      simulated_flows <- rtruncnorm(
        n    = length(fitted_values),
        a    = lower_bound,
        b    = upper_bound,
        mean = fitted_values,
        sd   = sd_vec
      )
      
      # 2) Matrix of simulated flows using dyads
      sim_df <- pairs_df
      sim_df$value <- simulated_flows
      flow_matrix <- create_flow_matrix(
        flows_df      = sim_df,
        country_names = country_names,
        value_col     = "value"
      )
      
      # 3) Matrix of observed flows for scaling
      obs_df <- pairs_df
      obs_df$value <- total_observed_flows
      obs_matrix <- create_flow_matrix(
        flows_df      = obs_df,
        country_names = country_names,
        value_col     = "value"
      )
      
      sim_row_sums <- rowSums(flow_matrix)
      obs_row_sums <- rowSums(obs_matrix)
      
      scaling_factors <- obs_row_sums / sim_row_sums
      
      # Row-wise scaling (R recycles across columns correctly here)
      flow_matrix <- flow_matrix * scaling_factors
      
      year_sims[[m]] <- flow_matrix
    }
    
    all_year_sims[[yr_chr]] <- year_sims
  }
  
  all_year_sims
}


# -------------------------------------------------------------------
# 4. Network measures
# -------------------------------------------------------------------

# Minimum reciprocity of a weighted network
reciprocity_min <- function(M){
  # Assumes empty diagonals
  lower <- M[lower.tri(M)]
  upper <- t(M)[lower.tri(M)]
  sum(pmin(lower, upper)) / sum(M)
}

# -------------------------------------------------------------------
# 5. Further helpers
# -------------------------------------------------------------------

compute_metrics <- function(sim_list) {
  tibble(
    reciprocity = map_dbl(sim_list, reciprocity_min),
    gini        = map_dbl(sim_list, migration.gini.total),
    acv         = map_dbl(sim_list, migration.acv),
    inequality  = map_dbl(
      sim_list,
      ~ migration.inequality(.x) / sum(rowSums(.x))
    )
  )
}


combine_metrics <- function(metrics_list, model_name) {
  bind_rows(metrics_list, .id = "Year") %>%
    group_by(Year) %>%
    mutate(sim_id = row_number()) %>%
    ungroup() %>%
    pivot_longer(
      cols      = c(reciprocity, gini, acv, inequality),
      names_to  = "statistic_raw",
      values_to = "value"
    ) %>%
    mutate(Model = model_name)
}
