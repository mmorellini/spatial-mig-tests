# SCRIPT DESCRIPTION -----------------------------------------------------------
# 01_load_clean_data.R
# Purpose: Load and harmonise bilateral migration flows, gravity covariates,
#          World Bank GDP/population, and UN migrant stocks into a single
#          panel dataset "merged" used in all PPML models and simulations.
# Inputs:
#   - data_raw/2024-10-15_HMigD_raw_results_table.xlsx
#   - data_raw/2023-07-15_CEPII-Gravity_V202211.rds
#   - World Bank API (GDP per capita PPP, population)
#   - UN migrant stocks via untools::getUNstocks(version = "2019")
# Output:
#   - data_processed/merged-flows.RData (object: merged)

library(dplyr)
library(tidyr)
library(readxl)
library(janitor)
library(stringr)
library(here)
library(countrycode)
library(wbstats)
library(readr)
library(untools)
library(mosaic)

# 1. Load bilateral migration data ---------------------------------------------

## Migration data comes from Maciej J Dańko, Arkadiusz Wiśniowski, Domantas
## Jasilionis, Dmitri A Jdanov, Emilio Zagheni, "Assessing the quality of data
## on international migration flows in Europe: The case of undercounting",
## Migration Studies, 12(2), 2024. https://doi.org/10.1093/migration/mnae014

## As of 15/10/2024, the data can be downloaded here:
## https://maciej-jan-danko.shinyapps.io/HMigD_Shiny_App_I/

flows <- read_excel(
  here("data_raw", "2024-10-15_HMigD_raw_results_table.xlsx")
)

# 2. Load and clean gravity data ----------------------------------------------

gravity <- readRDS(
  here::here("data_raw", "2023-07-15_CEPII-Gravity_V202211.rds")
)

## 2.1 Convert ISO3 to ISO2 to match HMigD coding
gravity <- gravity %>%
  mutate(
    sending_country   = countrycode(iso3_o,
                                    origin = "iso3c",
                                    destination = "iso2c"),
    receiving_country = countrycode(iso3_d,
                                    origin = "iso3c",
                                    destination = "iso2c")
  )

## Harmonise UK code
gravity$sending_country[gravity$sending_country == "GB"]   <- "UK"
gravity$receiving_country[gravity$receiving_country == "GB"] <- "UK"

## 2.2 Restrict to relevant years and add bilateral EU membership
gravity <- gravity %>%
  filter(year >= 2001, year < 2022) %>%
  mutate(
    eu_both = case_when(
      eu_o == 1 & eu_d == 1 ~ 1,
      .default = 0
    )
  )

## 2.3 Drop countries not present in HMigD flows and self-flows
gravity <- gravity %>%
  filter(
    sending_country %in% flows$sending_country,
    receiving_country %in% flows$receiving_country,
    sending_country != receiving_country
  )

## 2.4 Collapse any duplicates (e.g. due to German unification issues)
gravity <- gravity %>%
  group_by(sending_country, receiving_country, year) %>%
  fill(-c(sending_country, receiving_country, year), .direction = "updown") %>%
  slice(1) %>%
  ungroup()

# 3. Load World Bank GDP per capita and population -----------------------------

## Using World Bank API via wbstats
wb_gdp <- wb_data(
  country    = c(unique(flows$sending_country), "GB"),
  indicator  = "NY.GDP.PCAP.PP.CD",
  start_date = 2001,
  end_date   = 2021
) %>%
  clean_names() %>%
  rename(gdpcap_wb = ny_gdp_pcap_pp_cd)


## Plug in missing population data for Romania
wb_pop <- wb_data(
  country    = c(unique(flows$sending_country), "GB"),
  indicator  = "SP.POP.TOTL",
  start_date = 2001,
  end_date   = 2021
) %>%
  clean_names() %>%
  rename(pop_wb = sp_pop_totl)

## Merge GDP and population data and select variables of interest
wb_gdp <- wb_gdp %>%
  left_join(
    wb_pop,
    by = c("iso2c", "iso3c", "country", "date")
  ) %>%
  select(iso3c, country, date, gdpcap_wb, pop_wb)

# 4. Merge HMIgD flows, gravity, and World Bank data ---------------------------

## 4.1 Flows + gravity (join by country dyad and year)
merged <- flows %>%
  full_join(
    gravity,
    by = c("sending_country", "receiving_country", "year")
  )

## 4.2 Add GDP and population for origin and destination
merged <- merged %>%
  full_join(
    wb_gdp,
    by = c("year" = "date", "iso3_o" = "iso3c")
  ) %>%
  full_join(
    wb_gdp,
    by = c("year" = "date", "iso3_d" = "iso3c"),
    suffix = c("_o", "_d")
  ) %>%
  mutate(
    gdpcap_wb_o = gdpcap_wb_o,
    gdpcap_wb_d = gdpcap_wb_d,
    pop_wb_o    = pop_wb_o,
    pop_wb_d    = pop_wb_d
  )

# 5. Select variables of interest ----------------------------------------------

merged <- merged %>%
  select(
    sending_country, receiving_country,
    starts_with("pred"),
    pop_wb_o, pop_wb_d,
    dist, distcap, contig,
    sibling, sibling_ever, comcol, col45, col_dep_ever,
    comlang_off, comlang_ethno,
    eu_o, eu_d, eu_both,
    comleg_pretrans, comleg_posttrans,
    contains("_wb_"),
    year
  )

# 6. Add world regions and dyad identifiers ------------------------------------

## Use UN regions via countrycode::codelist
merged <- merged %>%
  left_join(
    countrycode::codelist %>%
      rename_all(\(x) paste0(x, "_o")) %>%
      mutate(
        sending_country = case_when(
          iso2c_o == "GB" ~ "UK",
          .default = iso2c_o
        )
      ) %>%
      select(
        sending_country, iso3c_o, iso.name.en_o,
        region_o, region23_o
      ),
    by = "sending_country"
  ) %>%
  left_join(
    countrycode::codelist %>%
      rename_all(\(x) paste0(x, "_d")) %>%
      mutate(
        receiving_country = case_when(
          iso2c_d == "GB" ~ "UK",
          .default = iso2c_d
        )
      ) %>%
      select(
        receiving_country, iso3c_d, iso.name.en_d,
        region_d, region23_d
      ),
    by = "receiving_country"
  ) %>%
  mutate(
    same_region = case_when(
      region23_o == region23_d ~ 1,
      .default = 0
    ),
    dyad = paste0(iso3c_o, iso3c_d)
  )

# 7. Add bilateral migrant stocks ----------------------------------------------

stocks <- untools::getUNstocks(version = "2019") %>%
  as_tibble() %>%
  filter(
    host_iso3   %in% merged$iso3c_d,
    origin_iso3 %in% merged$iso3c_o,
    year == 2000
  ) %>%
  mutate(dyad = paste0(origin_iso3, host_iso3))

stocks <- stocks %>% 
  rename(iso3c_o = origin_iso3,
         iso3c_d = host_iso3) %>% 
  select(year, iso3c_o, iso3c_d, dyad, stock)

## Expand all possible dyads x years 
stocks <- expand_grid(
  year    = sort(unique(stocks$year)),
  iso3c_o = sort(unique(stocks$iso3c_o)),
  iso3c_d = sort(unique(stocks$iso3c_d)) 
) %>%
  mutate(dyad = paste0(iso3c_o, iso3c_d)) %>%
  left_join(
    stocks %>% select(year, iso3c_o, iso3c_d, stock),
    by = c("year", "iso3c_o", "iso3c_d")
  ) %>%
  mutate(
    stock = tidyr::replace_na(stock, 0)
  ) %>% 
  filter(iso3c_o != iso3c_d) %>% 
  select(-year)

## Merge with flow data
merged <- stocks %>%
  full_join(merged)

# 8. Generate ratios and lagged variables --------------------------------------

merged <- merged %>%
  mutate(corridor = paste0(sending_country, "→", receiving_country),
         gdp_ratio = gdpcap_wb_o / gdpcap_wb_d,
         pop_ratio = pop_wb_d / pop_wb_o) %>%
  group_by(corridor) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    gdpcap_wb_o_lag = lag(gdpcap_wb_o, 1),
    gdpcap_wb_d_lag = lag(gdpcap_wb_d, 1),
    pop_wb_o_lag    = lag(pop_wb_o, 1),
    pop_wb_d_lag    = lag(pop_wb_d, 1),
    gdp_ratio_lag   = lag(gdp_ratio, 1),
    pop_ratio_lag   = lag(pop_ratio, 1),
    stocks_lag      = lag(stock, 1),
    eu_both_lag     = lag(eu_both, 1)
  ) %>%
  ungroup()

# Compute similarity indices

pop_rng_lag <- range(c(merged$pop_wb_o_lag, merged$pop_wb_d_lag),
                     na.rm = TRUE)
gdp_rng_lag <- range(c(merged$gdpcap_wb_o_lag, merged$gdpcap_wb_d_lag),
                     na.rm = TRUE)

merged <- merged %>%
  mutate(
    pop_sim_lag = 1 - abs(pop_wb_o_lag - pop_wb_d_lag) / (pop_rng_lag[2] - pop_rng_lag[1]),
    gdp_sim_lag = 1 - abs(gdpcap_wb_o_lag - gdpcap_wb_d_lag) / (gdp_rng_lag[2] - gdp_rng_lag[1])
  ) %>%
  mutate(
    # clamp for numerical safety
    pop_sim_lag   = pmax(pmin(pop_sim_lag, 1), 0),
    gdp_sim_lag   = pmax(pmin(gdp_sim_lag, 1), 0)
  )

# 9. Subset to period of interest ----------------------------------------------

merged <- merged %>% 
  filter(year >= 2002 & year < 2022)
  

# 10. Save merged file ---------------------------------------------------------

if (!dir.exists(here("data_processed"))) {
  dir.create(here("data_processed"), recursive = TRUE)
}

save(merged, file = here("data_processed", "merged-flows.RData"))
