# compute_mahalanobis.R
# Purpose:
#   - Compute Mahalanobis distance between observed and simulated flows
#     per country, per year, and per model (Gravity / MST).
#   - Compute Mahalanobis & other distance measures by corridor.
#
# Inputs:
#   - data_processed/obs_matrices_year.rds
#   - data_processed/gravity_year_sims.rds
#   - data_processed/mst_year_sims.rds

library(here)
library(dplyr)
library(purrr)
library(tidyr)

# ------------------------------------------------------------------------------
# Mahalanobis by country & year
# direction = "in", "out", or "both"
#   * Here we explicitly DROP the diagonal element (ctry -> ctry)
#     from both observed and simulated vectors before computing Mahalanobis.
# ------------------------------------------------------------------------------

get_mahalanobis_by_country <- function(obs, sim, model_name,
                                       tol = 1e-20, verbose = TRUE) {
  
  # Convert names to character to avoid mismatch issues
  obs_years <- as.character(names(obs))
  sim_years <- as.character(names(sim))
  
  # Find common years
  year_labels <- intersect(obs_years, sim_years)
  
  # Check if any years match
  if (length(year_labels) == 0) {
    stop("No matching years found between observed and simulated data. Check the names of your input lists.")
  }
  
  # Initialize list to store results
  result_list <- list()
  
  # Loop over years
  for (year in year_labels) {
    if (verbose) message("Processing year: ", year)
    
    # Extract observed migration matrix for the current year
    obs_matrix <- obs[[year]]
    
    # Extract simulated matrices for the current year
    sim_matrices <- sim[[year]]  # List of matrices
    
    # Get country names (column names = destinations)
    countries <- colnames(obs_matrix)
    
    # Initialize a list to store results for this year
    year_results <- list()
    
    # Loop over each country (destination)
    for (country in countries) {
      
      # 1. Observed migration flows for the country (column = inflows)
      obs_values <- obs_matrix[, country, drop = FALSE]
      
      # 2. All simulated flows for that country across simulations
      sim_values <- do.call(
        rbind,
        lapply(sim_matrices, function(mat) mat[, country, drop = FALSE])
      )
      
      # Ensure sim_values is a matrix
      sim_values <- as.matrix(sim_values)
      
      # 3. Mean and covariance matrix of simulated values
      center            <- colMeans(sim_values, na.rm = TRUE)
      covariance_matrix <- cov(sim_values, use = "complete.obs")
      
      # 4. Handle singular covariance matrices (minimal regularisation)
      if (det(covariance_matrix) < tol) {
        if (verbose) {
          message("Regularizing covariance matrix for year ", year,
                  " country ", country)
        }
        covariance_matrix <- covariance_matrix + diag(tol, nrow(covariance_matrix))
      }
      
      # 5. Mahalanobis distance
      mhd <- mahalanobis(
        x   = obs_values,
        center = center,
        cov    = covariance_matrix,
        tol    = tol
      )
      
      # Store in results list
      year_results[[country]] <- data.frame(
        year        = as.integer(year),
        Year        = as.integer(year),  # convenient numeric column
        country     = country,
        model       = model_name,
        mahalanobis = as.numeric(mhd),
        stringsAsFactors = FALSE
      )
    }
    
    # Combine all results for this year
    result_list[[year]] <- do.call(rbind, year_results)
  }
  
  # Combine all years into a single dataframe & order
  mahalanobis_df <- do.call(rbind, result_list)
  mahalanobis_df <- mahalanobis_df[order(mahalanobis_df$year,
                                         mahalanobis_df$country), ]
  
  rownames(mahalanobis_df) <- NULL
  return(mahalanobis_df)
}

# ------------------------------------------------------------------------------
# Mahalanobis by corridor (origin-destination dyad)
#   corridors: data.frame/tibble with columns origin, destination (ISO2 codes)
#   if corridors = NULL, we do all dyads (except i = j)
#   Note: here we already exclude diagonal by origin != destination.
# ------------------------------------------------------------------------------

get_mahalanobis_by_corridor <- function(obs, sim, model_name,
                                        tol = 1e-20, verbose = TRUE,
                                        corridors = NULL) {
  
  # Years as character
  obs_years <- as.character(names(obs))
  sim_years <- as.character(names(sim))
  year_labels <- intersect(obs_years, sim_years)
  
  if (length(year_labels) == 0) {
    stop("No matching years found between observed and simulated data. Check the names of your input lists.")
  }
  
  year_labels <- sort(year_labels)
  
  # If corridors not supplied, use all ordered OD pairs with i != j
  if (is.null(corridors)) {
    sample_mat    <- obs[[year_labels[1]]]
    all_countries <- colnames(sample_mat)
    
    corridors <- expand.grid(
      origin      = all_countries,
      destination = all_countries,
      stringsAsFactors = FALSE
    )
    corridors <- corridors[corridors$origin != corridors$destination, ]
  }
  
  result_list <- list()
  
  # Loop over years
  for (year in year_labels) {
    if (verbose) message("Processing year (corridors): ", year)
    
    obs_matrix  <- obs[[year]]
    sim_matrices <- sim[[year]]  # list of matrices for this year
    
    year_rows <- vector("list", nrow(corridors))
    
    # Loop over corridors
    for (i in seq_len(nrow(corridors))) {
      o <- corridors$origin[i]
      d <- corridors$destination[i]
      
      # Observed value for corridor o -> d
      x_obs <- obs_matrix[o, d]
      
      # Simulated values for that corridor in this year
      x_sim <- vapply(
        sim_matrices,
        function(mat) mat[o, d],
        numeric(1)
      )
      
      mu <- mean(x_sim, na.rm = TRUE)
      v  <- var(x_sim, na.rm = TRUE)
      
      # 1D Mahalanobis (squared z-score); guard against zero variance
      if (is.na(v) || v < tol) {
        mhd <- 0
      } else {
        mhd <- (x_obs - mu)^2 / v
      }
      
      year_rows[[i]] <- data.frame(
        year        = as.integer(year),
        Year        = as.integer(year),
        origin      = o,
        destination = d,
        model       = model_name,
        mahalanobis = as.numeric(mhd),
        stringsAsFactors = FALSE
      )
    }
    
    result_list[[year]] <- do.call(rbind, year_rows)
  }
  
  mahalanobis_corr <- do.call(rbind, result_list)
  mahalanobis_corr <- mahalanobis_corr[order(mahalanobis_corr$year,
                                             mahalanobis_corr$origin,
                                             mahalanobis_corr$destination), ]
  rownames(mahalanobis_corr) <- NULL
  
  return(mahalanobis_corr)
}
