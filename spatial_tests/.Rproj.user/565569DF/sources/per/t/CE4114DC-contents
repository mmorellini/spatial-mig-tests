# SCRIPT DESCRIPTION -----------------------------------------------------------
# 09_ppml_eurostat_data.R
# Purpose: Fit PPML models with Eurostat migration flow data as outcome:
#   - Baseline (FE only)
#   - Gravity (distance only)
#   - MST (full set of covariates with lags)
# Perform RESET-style tests and export LaTeX regression table.
# Inputs:
#   - data_processed/merged-flows.RData (merged); for country codes & covariates
#   - Eurostat migration flow data (MIGR_IMM5PRV), accessed via API
# Output:
#   - tables/ppml_eurostat_data.tex

library(dplyr)
library(tidyr)
library(fixest)
library(here)
library(eurostat)

# 0. Load merged data and Eurostat migration flow data -------------------------

source(here("scripts", "utils", "plotting_helpers.R"))

if (!exists("merged")) {
  load(here("data_processed", "merged-flows.RData"))
}

## define country codes of chosen countries
chosen_country_codes <- unique(c(merged$sending_country,
                                 merged$receiving_country))

## Get Eurostat data with eurostat package and data code MIGR_IMM5PRV
eurostat_df <- eurostat::get_eurostat("migr_imm5prv", time_format = "num")

eurostat_df$geo[eurostat_df$geo == "EL"] <- "GR"
eurostat_df$partner[eurostat_df$partner == "EL"] <- "GR"

eurostat_df <- eurostat_df %>% 
  filter(TIME_PERIOD > 2001 & TIME_PERIOD <= 2021 & 
           age == "TOTAL" &
           agedef == "REACH" &
           sex == "T" &
           geo %in% chosen_country_codes &
           partner %in% chosen_country_codes
           ) %>% 
  filter(partner %in% unique(eurostat_df$geo) &
           partner != geo) %>% 
  rename(year = TIME_PERIOD,
         sending_country = partner,
         receiving_country = geo)

## Add covariates from merged data
eurostat_df <- eurostat_df %>%
  left_join(
    merged %>%
      select(
        sending_country,
        receiving_country,
        year,
        distcap,
        pop_ratio_lag,
        contig,
        comlang_off,
        gdp_ratio_lag,
        stocks_lag,
        sibling_ever,
        eu_both_lag
      ) %>%
      distinct(),
    by = c("sending_country", "receiving_country", "year")
  ) %>% 
  drop_na()


# 1. Baseline model ------------------------------------------------------------

base_pois <- fepois(
  values ~ 1 |
    (sending_country + receiving_country) + year,
  vcov = "threeway",
  data = eurostat_df
)

# 2. Gravity model -------------------------------------------------------------

gravity_pois <- fepois(
  values ~ log1p(distcap) + log1p(pop_ratio_lag) |
    (sending_country + receiving_country) + year,
  vcov = "threeway",
  data = eurostat_df
)

# 3. MST model -----------------------------------------------------------------

mst_pois <- fepois(
  values ~ log1p(distcap) +
    log1p(pop_ratio_lag) +
    contig +
    comlang_off +
    log1p(gdp_ratio_lag) +
    log1p(stocks_lag) +
    sibling_ever +
    eu_both_lag |
    (sending_country + receiving_country) + year,
  data = eurostat_df,
  vcov = "threeway"
)

# 4. RESET-style tests ---------------------------------------------------------

## Construct squared linear predictors for each model
eurostat_df <- eurostat_df %>%
  mutate(
    predict_base    = (predict(base_pois,    type = "link"))^2,
    predict_gravity = (predict(gravity_pois, type = "link"))^2,
    predict_mst     = (predict(mst_pois,     type = "link"))^2
  )

## Baseline RESET
base_reset <- fepois(
  values ~ 1 + predict_base |
    (sending_country + receiving_country) + year,
  data = eurostat_df,
  vcov = "threeway"
)

## Gravity RESET
gravity_reset <- fepois(
  values ~ log1p(distcap) + log1p(pop_ratio_lag) + predict_gravity |
    (sending_country + receiving_country) + year,
  data = eurostat_df,
  vcov = "threeway"
)

## MST RESET
mst_reset <- fepois(
  values ~ log1p(distcap) +
    log1p(pop_ratio_lag) +
    contig +
    comlang_off +
    log1p(gdp_ratio_lag) +
    log1p(stocks_lag) +
    sibling_ever +
    eu_both_lag +
    predict_mst |
    (sending_country + receiving_country) + year,
  data = eurostat_df,
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
  file = here("tables", "ppml_eurostat_data.tex"),
  extralines = list(
    "__Ramsey RESET test (p-value)" = res
  )
)

# 6. Eurostat inflows and outflows plot ----------------------------------------

flows_summary <- eurostat_df %>%
  group_by(sending_country, year) %>%
  summarise(
    outflows = sum(values, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    eurostat_df %>%
      group_by(receiving_country, year) %>%
      summarise(
        inflows = sum(values, na.rm = TRUE),
        .groups     = "drop"
      ),
    by = c("sending_country" = "receiving_country", "year")
  ) %>%
  rename(Country = sending_country) %>%
  mutate(
    Country = if_else(Country == "UK", "GB", Country),
    Country = countrycode(
      Country,
      origin      = "iso2c",
      destination = "country.name"
    )
  ) %>%
  pivot_longer(
    cols      = -c(Country, year),
    names_to  = "Flow",
    values_to = "Value"
  ) %>% 
  mutate(Year     = as.numeric(year),
         Value    = as.numeric(Value),
         Flow = case_when(Flow == "inflows"  ~ "Inflows",
                          Flow == "outflows" ~ "Outflows",
                          TRUE               ~ Flow
    ),
    Flow = factor(Flow, levels = c("Inflows", "Outflows"))
  )


flow_colors <- c(
  "Inflows"  = "#CC79A7",  # Reddish purple
  "Outflows" = "#009E73"   # Bluish green
)

flow_linetypes <- c(
  "Inflows"  = "solid",
  "Outflows" = "dashed"
)

p_migration_flows <- ggplot(
  flows_summary,
  aes(
    x        = Year,
    y        = Value / 1000,
    color    = Flow,
    fill     = Flow,
    linetype = Flow
  )
) +
  facet_wrap(
    ~ Country,
    scales = "free_y",
    ncol   = 4,
    strip.position = "top"
  ) +
  # Blended ribbons for IQR
   geom_line(linewidth = 1.2) +
  scale_color_manual(values = flow_colors) +
  scale_fill_manual(values  = flow_colors) +
  scale_linetype_manual(values = flow_linetypes) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 4),
    labels = label_number(accuracy = 1)
  ) +
  theme_ready_pub(base_size = 18) +
  labs(
    title = "Eurostat Migration Inflows & Outflows in Europe by Country, 2002–2021",
    subtitle = "Y-axis scales vary by country.",
    x = "Year",
    y = "Migration volume (thousands)",
    caption = "Source: Eurostat (2023), Immigration by Age Group, Sex and Country of Previous Residence (MIGR_IMM5PRV).",
    color   = "Flow type:",
    fill    = "Flow type:",
    linetype = "Flow type:"
  ) +
  theme(
    legend.position  = "right",
    legend.direction = "vertical",
    legend.key.width  = unit(2.5, "cm"),
    legend.key.height = unit(0.7, "cm")
  )

p_migration_flows <- shift_legend2(p_migration_flows)

ggsave(
  p_migration_flows,
  file   = here("figures", "migration_flows_eurostat.png"),
  device = png,
  type   = "cairo",
  width  = 15,
  height = 18,
  units = "in",
  dpi    = 300
)
