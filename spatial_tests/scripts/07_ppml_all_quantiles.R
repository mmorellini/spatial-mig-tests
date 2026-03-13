# SCRIPT DESCRIPTION -----------------------------------------------------------
# 07_ppml_all_quantiles.R
# Purpose: Fit PPML gravity and MST models for all quantile outcomes:
#   - pred_q05, pred_q25, pred_q50, pred_q75, pred_q95
# Specification matches main q50 models.
# Inputs:
#   - data_processed/merged-flows.RData (object: merged)
# Output:
#   - tables/ppml_gravity_all_quantiles.tex
#   - tables/ppml_mst_all_quantiles.tex

library(dplyr)
library(fixest)
library(here)

# 0. Load processed data -------------------------------------------------------

if (!exists("merged")) {
  load(here("data_processed", "merged-flows.RData"))
}

# 1. Define outcomes and names -------------------------------------------------

outcomes <- c("pred_q05", "pred_q25", "pred_q50", "pred_q75", "pred_q95")
qnames   <- c("q05",      "q25",      "q50",      "q75",      "q95")


# 2. Baseline models for all quantiles -----------------------------------------

base_pois_all <- lapply(outcomes, function(y) {
  fepois(
    as.formula(
      paste0(
        y,
        " ~ 1 | (sending_country + receiving_country) + year"
      )
    ),
    vcov = "threeway",
    data = merged
  )
})

names(base_pois_all) <- qnames

# 3. Gravity models for all quantiles -----------------------------------------

gravity_pois_all <- lapply(outcomes, function(y) {
  fepois(
    as.formula(
      paste0(
        y,
        " ~ log1p(distcap) + log1p(pop_ratio_lag) | (sending_country + receiving_country) + year"
      )
    ),
    vcov = "threeway",
    data = merged
  )
})

names(gravity_pois_all) <- qnames

# 4. MST models for all quantiles ---------------------------------------------

mst_pois_all <- lapply(outcomes, function(y) {
  fepois(
    as.formula(
      paste0(
        y,
        " ~ log1p(distcap) + log1p(pop_ratio_lag) + contig + comlang_off + ",
        "log1p(gdp_ratio_lag) + log1p(stocks_lag) + ",
        "sibling_ever + eu_both_lag | ",
        "(sending_country + receiving_country) + year"
      )
    ),
    vcov = "threeway",
    data = merged
  )
})

names(mst_pois_all) <- qnames

# 5. Ensure tables directory exists -------------------------------------------

if (!dir.exists(here("tables"))) {
  dir.create(here("tables"), recursive = TRUE)
}

# 5. Export LaTeX tables -------------------------------------------------------


# 5.2 Baseline models
etable(
  base_pois_all,
  fitstat    = c("n", "cor2", "apr2", "wapr2"),
  vcov       = "threeway",
  signif.code = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "." = 0.10),
  file       = here("tables", "ppml_base_all_quantiles.tex"),
  tex        = TRUE,
  replace = TRUE
)

# 5.2 Gravity models
etable(
  gravity_pois_all,
  fitstat    = c("n", "cor2", "apr2", "wapr2"),
  vcov       = "threeway",
  signif.code = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "." = 0.10),
  file       = here("tables", "ppml_gravity_all_quantiles.tex"),
  tex        = TRUE,
  replace = TRUE
)

# 5.3 MST models
etable(
  mst_pois_all,
  fitstat    = c("n", "cor2", "apr2", "wapr2"),
  vcov       = "threeway",
  signif.code = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "." = 0.10),
  file       = here("tables", "ppml_mst_all_quantiles.tex"),
  tex        = TRUE,
  replace = TRUE
)
