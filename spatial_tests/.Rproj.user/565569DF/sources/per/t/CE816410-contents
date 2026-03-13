# SCRIPT DESCRIPTION -----------------------------------------------------------
# 08_ppml_exclude_covid_q50.R
# Purpose: Fit PPML models with pred_q50 as outcome, excluding 2019-2021 
# (years of COVID-19 pandemic):
#   - Baseline (FE only)
#   - Gravity (distance only)
#   - MST (full set of covariates with lags)
# Perform RESET-style tests and export LaTeX regression table.
# Inputs:
#   - data_processed/merged-flows.RData (object: merged), filtered to years < 2019
# Output:
#   - tables/ppml_exclude_covid_q50.tex

library(dplyr)
library(tidyr)
library(fixest)
library(here)

# 0. Load processed data and filter years --------------------------------------

if (!exists("merged")) {
  load(here("data_processed", "merged-flows.RData"))
}

merged <- merged %>% 
  filter(year < 2019)

# 1. Baseline model ------------------------------------------------------------

base_pois <- fepois(
  pred_q50 ~ 1 |
    (sending_country + receiving_country) + year,
  vcov = "threeway",
  data = merged
)

# 2. Gravity model -------------------------------------------------------------

gravity_pois <- fepois(
  pred_q50 ~ log1p(distcap) + log1p(pop_ratio_lag) |
    (sending_country + receiving_country) + year,
  vcov = "threeway",
  data = merged
)

# 3. MST model -----------------------------------------------------------------

mst_pois <- fepois(
  pred_q50 ~ log1p(distcap) +
    log1p(pop_ratio_lag) +
    contig +
    comlang_off +
    log1p(gdp_ratio_lag) +
    log1p(stocks_lag) +
    sibling_ever +
    eu_both_lag |
    (sending_country + receiving_country) + year,
  data = merged,
  vcov = "threeway"
)

# 4. RESET-style tests ---------------------------------------------------------

## Construct squared linear predictors for each model
merged <- merged %>%
  mutate(
    predict_base    = (predict(base_pois,    type = "link"))^2,
    predict_gravity = (predict(gravity_pois, type = "link"))^2,
    predict_mst     = (predict(mst_pois,     type = "link"))^2
  )

## Baseline RESET
base_reset <- fepois(
  pred_q50 ~ 1 + predict_base |
    (sending_country + receiving_country) + year,
  data = merged,
  vcov = "threeway"
)

## Gravity RESET
gravity_reset <- fepois(
  pred_q50 ~ log1p(distcap) + log1p(pop_ratio_lag) + predict_gravity |
    (sending_country + receiving_country) + year,
  data = merged,
  vcov = "threeway"
)

## MST RESET
mst_reset <- fepois(
  pred_q50 ~ log1p(distcap) +
    log1p(pop_ratio_lag) +
    contig +
    comlang_off +
    log1p(gdp_ratio_lag) +
    log1p(stocks_lag) +
    sibling_ever +
    eu_both_lag +
    predict_mst |
    (sending_country + receiving_country) + year,
  data = merged,
  vcov = "threeway"
)

## Extract RESET p-values (for squared prediction term)
res <- c(
  base_reset$coeftable["predict_base",    "Pr(>|z|)"],
  gravity_reset$coeftable["predict_gravity", "Pr(>|z|)"],
  mst_reset$coeftable["predict_mst",      "Pr(>|z|)"])

# 5. Export LaTeX table --------------------------------------------------------

if (!dir.exists(here("tables"))) {
  dir.create(here("tables"), recursive = TRUE)
}

etable(
  base_pois, gravity_pois, mst_pois,
  fitstat = c("n", "cor2", "apr2", "wapr2"),
  vcov = "threeway",
  signif.code = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "." = 0.10),
  tex  = TRUE,
  replace = TRUE,
  file = here("tables", "ppml_exclude_covid_q50.tex"),
  extralines = list(
    "__Ramsey RESET test (p-value)" = res
  )
)
