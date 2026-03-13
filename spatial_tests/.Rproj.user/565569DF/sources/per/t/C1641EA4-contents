# SCRIPT DESCRIPTION -----------------------------------------------------------
# scripts/06_make_plots.R
# Purpose:
#   - Generate all plots contained the paper:
#       1) Migration inflows and outflows by country over time
#       2) Global network metrics over time (ACV, Gini, inequality, reciprocity)
#       3) Mahalanobis distances of inflows by country and over time
#       4) Scatterplot of simulated vs observed corridor flows, with outlier labels
#       5) Ribbon plots of observed vs simulated flows over tome in selected corridors
# Inputs:
#   - data_processed/merged-flows.RData
#   - data_processed/gravity_year_sims.rds
#   - data_processed/mst_year_sims.rds
#   - data_processed/global_values_long.rds
#   - data_processed/observed_values_all_q.rds
#   - data_processed/country_distances_long.rds
#   - data_processed/sim_summaries.rds
# Outputs:
#   - figures/migration_flows_hmigd.png
#   - figures/migration_metrics_trend.png
#   - figures/mahalanobis_country_time.png
#   - figures/scatter_sim_vs_obs.png
#   - figures/selected_corridors_flows.png


library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(countrycode)
library(scales)
library(readr)
library(purrr)
library(ggblend)
library(ggrepel)
library(viridis)


# 0. Load script with plotting helpers and data objects-------------------------

## Plotting helpers

source(here("scripts", "utils", "plotting_helpers.R"))

# Migration flows (Bayesian)
if (!exists("merged")) {
  load(here("data_processed", "merged-flows.RData"))  # loads `merged`
}

# Simulations from script 04
if (!exists("gravity_year_sims")) {
gravity_year_sims <- readRDS(
  here("data_processed", "gravity_year_sims.rds"))
}

if (!exists("mst_year_sims")) {
mst_year_sims <- readRDS(
  here("data_processed", "mst_year_sims.rds"))
}

# Network metrics from script 04
if (!exists("global_values_long")) {
global_values_long <- readRDS(
  here("data_processed", "global_values_long.rds"))
}

# Load observed indices for all pred_q* from script 04
if (!exists("obs_values_all_q")) {
obs_values_all_q <- readRDS(
  here("data_processed", "observed_values_all_q.rds"))
}

# By-country distances from script 05
if (!exists("country_distances_long")) {
country_distances_long <- readRDS(
  here("data_processed", "country_distances_long.rds"))
}

# By corridor summaries from script 05
if (!exists("sim_summaries")) {
sim_summaries <- readRDS(
  here("data_processed", "sim_summaries.rds"))
}

if (!exists("corridor_fit")) {
  corridor_fit <- readRDS(
    here("data_processed", "sim_summaries.rds"))
}


# 1. Migration inflows and outflows, Bayesian ribbons, faceted by country-------

flows_summary <- merged %>%
  group_by(sending_country, year) %>%
  summarise(
    q25_outflows = sum(pred_q25, na.rm = TRUE),
    q50_outflows = sum(pred_q50, na.rm = TRUE),
    q75_outflows = sum(pred_q75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    merged %>%
      group_by(receiving_country, year) %>%
      summarise(
        q25_inflows = sum(pred_q25, na.rm = TRUE),
        q50_inflows = sum(pred_q50, na.rm = TRUE),
        q75_inflows = sum(pred_q75, na.rm = TRUE),
        .groups     = "drop"
      ),
    by = c("sending_country" = "receiving_country", "year")
  ) %>%
  # Long: inflows/outflows + quantiles
  pivot_longer(
    cols      = -c(sending_country, year),
    names_to  = c("Quantile", "Flow"),
    names_pattern = "q(\\d+)_(.*)",
    values_to = "Value"
  ) %>%
  mutate(
    Year     = as.numeric(year),
    Value    = as.numeric(Value),
    Quantile = as.numeric(Quantile)
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
  # Wide again: q25 / q50 / q75 columns
  pivot_wider(
    names_from  = Quantile,
    values_from = Value,
    names_prefix = "q"
  ) %>%
  mutate(
    Flow = case_when(
      Flow == "inflows"  ~ "Inflows",
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
  (
    geom_ribbon(
      aes(
        ymin = q25 / 1000,
        ymax = q75 / 1000
      ),
      alpha     = 0.3,
      linewidth = 0.4
    ) %>% 
      blend("multiply")
  ) +
  # Median lines
  geom_line(aes(y = q50 / 1000), linewidth = 1.2) +
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
    title = "Estimated Migration Inflows & Outflows in Europe by Country, 2002–2021",
    subtitle = "Bayesian estimates showing median values (lines) and interquartile ranges (ribbons).\nY-axis scales vary by country.",
    x = "Year",
    y = "Migration volume (thousands)",
    caption = "Source: Human Migration Database (2023).",
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
  file   = here("figures", "migration_flows_hmigd.png"),
  device = png,
  type   = "cairo",
  width  = 15,
  height = 18,
  units = "in",
  dpi    = 300
)

# 2. Global network metrics over time (ACV, Gini, inequality, reciprocity)------

global_values_long <- global_values_long %>%
  mutate(Year = as.numeric(Year))

global_values_long_summary <- global_values_long %>%
  group_by(Year, statistic_raw, statistic, Model) %>%
  summarise(
    mean_value = median(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    lower      = min(value, na.rm = TRUE),
    upper      = max(value, na.rm = TRUE),
    .groups    = "drop"
  )

observed_values <- global_values_long %>%
  select(Year, statistic_raw, statistic, obs_value) %>%
  distinct() %>%
  mutate(
    Model      = "Observed",
    mean_value = obs_value
  ) %>%
  select(Year, statistic_raw, statistic, Model, mean_value)

global_values_long_summary <- bind_rows(
  global_values_long_summary,
  observed_values
) %>%
  mutate(
    statistic = factor(
      statistic,
      levels = c("ACV", "Gini", "Migration inequality", "Reciprocity")
    )
  )

# Keep only q25, q50, q75 and reshape to wide to get IQR
obs_iqr <- obs_values_all_q %>%
  filter(quantile %in% c("q25", "q50", "q75")) %>%
  select(Year, statistic_raw, quantile, obs_value) %>%
  mutate(Year = as.numeric(Year)) %>%
  tidyr::pivot_wider(
    names_from  = quantile,
    values_from = obs_value
  ) %>%
  rename(
    obs_q25 = q25,
    obs_q50 = q50,
    obs_q75 = q75
  ) %>%
  mutate(
    Model = "Observed",
    statistic = recode(
      statistic_raw,
      acv         = "ACV",
      gini        = "Gini",
      inequality  = "Migration inequality",
      reciprocity = "Reciprocity"
    ),
    statistic = factor(
      statistic,
      levels = c("ACV", "Gini", "Migration inequality", "Reciprocity")
    )
  )

# Colors and linetypes
model_colors <- c(
  "Observed" = "#000000",
  "Gravity"  = "#0072B2",
  "MST"      = "#E69F00"
)
model_fills <- c(
  "Gravity"  = "#A6CEE3",
  "MST"      = "#FDBF6F",
  "Observed" = "#3B3B3B"
)
model_linetypes <- c(
  "Observed" = "solid",
  "Gravity"  = "dotted",
  "MST"      = "longdash"
)

# Ensure factors are consistently defined
global_values_long_summary <- global_values_long_summary %>%
  mutate(
    Model = factor(Model, levels = c("Observed", "Gravity", "MST")),
    statistic = factor(
      statistic,
      levels = c("ACV", "Gini", "Migration inequality", "Reciprocity")
    )
  )

obs_iqr <- obs_iqr %>%
  mutate(
    Model = factor("Observed", levels = c("Observed", "Gravity", "MST"))
  )

# Ribbons for observed (IQR)
obs_ribbons <- obs_iqr %>%
  transmute(
    Year,
    statistic,
    Model,
    ymin = obs_q25,
    ymax = obs_q75
  )

# Ribbons for simulated models (Gravity & MST)
model_ribbons <- global_values_long_summary %>%
  filter(Model != "Observed") %>%
  transmute(
    Year,
    statistic,
    Model,
    ymin = lower,
    ymax = upper
  )

# Unified ribbons df for blending
ribbons_df <- bind_rows(obs_ribbons, model_ribbons)

p_metrics_trend <- ggplot(
  global_values_long_summary,
  aes(
    x        = Year,
    y        = mean_value,
    color    = Model,
    fill     = Model,
    linetype = Model
  )
) +
  facet_wrap(
    ~ statistic,
    scales = "free_y",
    ncol   = 2,
    strip.position = "top"
  ) +
  # One blended ribbon layer for Observed + Gravity + MST
  (
    geom_ribbon(
      data = ribbons_df,
      aes(
        x     = Year,
        ymin  = ymin,
        ymax  = ymax,
        fill  = Model,
        colour = Model,
        linetype = Model,
        group = interaction(Model, statistic)
      ),
      inherit.aes = FALSE,
      alpha       = 0.4,
      linewidth   = 0.4
    ) %>% 
      partition(vars(Model)) %>%    # blend per model
      blend("multiply")
  ) +
  # Lines for observed + simulated means
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values  = model_fills) +
  scale_linetype_manual(values = model_linetypes) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 4),
    labels = number_format(accuracy = 0.01)
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  guides(
    color    = guide_legend(nrow = 1),
    fill     = guide_legend(nrow  = 1),
    linetype = guide_legend(
    nrow     = 1
    )
  ) +
  theme_ready_pub(base_size = 18) +
  labs(
    title = "Observed and Simulated Migration Indices, 2002–2021",
    subtitle = paste0(
      "Observed values in black, with median values (lines) and interquartile ranges (ribbons).\n",
      "Simulated values in blue (gravity) and orange (MST), with simulation means (lines)\n",
      "and entire range of distributions (ribbons). Y-axis scales differ by index."
    ),
    x       = "Year",
    y       = "Index value",
    color   = "Model:",
    fill    = "Model:",
    linetype = "Model:"
  ) +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.key.width  = unit(2.5, "cm"),
    legend.key.height = unit(0.7, "cm")
  )

ggsave(
  p_metrics_trend,
  file   = here("figures", "migration_metrics_trend.png"),
  device = png,
  type   = "cairo",
  width  = 13,
  height = 13,
  units = "in",
  dpi    = 300
)

# 3. Mahalanobis distances over time, faceted by country------------------------
 
mahalanobis_df <- country_distances_long %>%
  mutate(
    Year  = as.numeric(Year),
    Model = case_when(
      model == "gravity" ~ "Gravity",
      model == "mst"     ~ "MST",
      TRUE               ~ model
    ),
    Country = if_else(country == "UK", "GB", country),
    Country = countrycode(
      Country,
      origin      = "iso2c",
      destination = "country.name"
    )
  )

mah_colors <- c("Gravity" = "#0072B2", "MST" = "#E69F00")
mah_linetypes <- c("Gravity" = "solid", "MST" = "solid")

p_mahalanobis <- ggplot(
  mahalanobis_df,
  aes(
    x        = Year,
    y        = mahalanobis,
    color    = Model,
    linetype = Model
  )
) +
  facet_wrap(
    ~ Country,
    scales = "free_y",
    ncol   = 4,
    strip.position = "top"
  ) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values    = mah_colors) +
  scale_linetype_manual(values = mah_linetypes) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  guides(
    color    = guide_legend(nrow = 2),
    linetype = guide_legend(
      #override.aes = list(linewidth = 1.3),
      nrow         = 2
    )
  ) +
  theme_ready_pub(base_size = 18) +
  labs(
    title = "Distance Between Observed and Simulated Inflows by Country, 2002–2021",
    subtitle = paste0("Mahalanobis distance between observed and simulated inflows, with gravity values in blue\n",
                      "and MST values in orange. Y-axis scales vary by country."),
    x       = "Year",
    y       = "Mahalanobis distance",
    color   = "Model:",
    linetype = "Model:"
  ) +
  theme(legend.position  = "right",
        legend.direction = "vertical",
        legend.key.width  = unit(2.5, "cm"),
        legend.key.height = unit(0.7, "cm"))

p_mahalanobis <- shift_legend2(p_mahalanobis)

ggsave(
  p_mahalanobis,
  file   = here("figures", "mahalanobis_country_time.png"),
  device = png,
  type   = "cairo",
  width  = 15,
  height = 18,
  units = "in",
  dpi    = 300
)


# 4. Scatterplot of observed vs simulated corridor-years observations-----------

# Generate labels for example outliers:
corridor_fit_extreme <- corridor_fit %>%
  mutate(origin_name = countrycode::countrycode(iso3c_o,
                                                "iso3c",
                                                "country.name"),
         dest_name   = countrycode::countrycode(iso3c_d,
                                                "iso3c",
                                                "country.name")) %>% 
  mutate(label = paste0(origin_name, " → ", dest_name, " \n(", year, ")")) %>% 
  filter(label %in% c("Estonia → Finland \n(2013)",
                      "Finland → Estonia \n(2017)",
                      "Cyprus → Greece \n(2021)",
                      "Austria → Germany \n(2004)"))

# Create unique IDs for highlight condition
corridor_fit <- corridor_fit %>%
  mutate(id = paste(origin, destination, year))

extreme_ids <- corridor_fit_extreme %>%
  mutate(id = paste(origin, destination, year)) %>%
  pull(id)

# Shapes for distinction
# Shapes for background points (grey)
bg_shapes <- c(
  "Gravity" = 21,  # hollow circle
  "MST"     = 24   # hollow triangle
)

p_sim_vs_obs <- ggplot(corridor_fit, aes(
  x = pred_q50,
  y = mean_sim
)) +
  # Dashed perfect-fit line
  geom_abline(
    slope     = 1,
    intercept = 0,
    linetype  = "dashed",
    colour    = "black"
  ) +
  
  # 1) All points, grey, with model-specific shapes 21 / 24
  geom_point(
    aes(shape = model),
    colour = "grey70",
    alpha  = 0.2,
    size   = 2
  ) +
  scale_shape_manual(
    values = bg_shapes,
    guide  = "none"   # don't show these shapes in the legend
  ) +
  
  # 2) Highlighted points: coloured, fully opaque, shapes 16 (Gravity) and 17 (MST)
  geom_point(
    data = corridor_fit %>% filter(id %in% extreme_ids, model == "Gravity"),
    aes(colour = model),
    shape = 16,   # filled circle
    alpha = 1,
    size  = 4
  ) +
  geom_point(
    data = corridor_fit %>% filter(id %in% extreme_ids, model == "MST"),
    aes(colour = model),
    shape = 17,   # filled triangle
    alpha = 1,
    size  = 4
  ) +
  
  # 3) Labels only on highlighted points
  geom_label_repel(
    data        = corridor_fit_extreme,
    inherit.aes = FALSE,
    aes(
      x     = pred_q50,
      y     = mean_sim,
      label = label,
      colour = model
    ),
    size              = 5,
    fontface          = "bold",
    box.padding       = 0.4,
    label.size        = 0.4,
    max.overlaps      = Inf,
    min.segment.length = 0.1,
    fill              = "white",
    show.legend       = FALSE
  ) +
  
  # Theme & scales
  theme_ready_pub(base_size = 18) +
  scale_colour_manual(values = c(
    "Gravity" = "#0072B2",
    "MST"     = "#E69F00"
  )) +
  scale_x_continuous(
    trans  = "log10",
    labels = label_number(accuracy = 1),
    breaks = 10^(0:5)
  ) +
  scale_y_continuous(
    trans  = "log10",
    labels = label_number(accuracy = 1),
    breaks = 10^(0:5)
  ) +
  
  labs(
    x      = "Observed flows (log scale)",
    y      = "Simulated flows (mean across simulations, log scale)",
    colour = "Model:",
    title  = "Simulated vs Observed Flows by Corridor–Year and Model",
    subtitle = paste0(
      "Corridor–year observations for gravity (circles) and MST (triangles) in grey.\n",
      "Selected corridors highlighted and coloured in blue (gravity) and orange (MST)."
    )
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 6, alpha = 1))
  )

ggsave(
  p_sim_vs_obs,
  file   = here("figures", "scatter_sim_vs_obs.png"),
  device = png,
  type   = "cairo",
  width  = 13,
  height = 13,
  units = "in",
  dpi    = 300
)

# 5. Ribbon plot: selected corridors (Observed vs Gravity vs MST)---------------

# Select corridors to plot from simulations
chosen_corridors <- sim_summaries %>% 
  mutate(origin_name = countrycode::countrycode(iso3c_o, "iso3c", "country.name"),
         dest_name   = countrycode::countrycode(iso3c_d, "iso3c", "country.name")) %>% 
  mutate(corridor_name = paste(origin_name, "→", dest_name)) %>% 
  mutate(corridor = factor(corridor_name)) %>% 
  filter(corridor_pair %in% c("CY → GR",
                              "EE → FI",
                              "FI → EE",
                              "AT → DE"))

# Observed (Bayesian) corridor series from `merged`
obs_corr <- merged %>%
  mutate(origin_name = countrycode::countrycode(iso3c_o, "iso3c", "country.name"),
         dest_name   = countrycode::countrycode(iso3c_d, "iso3c", "country.name")) %>% 
  mutate(corridor_name = paste(origin_name, "→", dest_name)) %>% 
  semi_join(
    chosen_corridors %>% select(origin, destination, corridor_name),
    by = c("sending_country" = "origin",
           "receiving_country" = "destination",
           "corridor_name" = "corridor_name")
  ) %>%
  transmute(
    Year          = as.numeric(year),
    origin        = sending_country,
    destination   = receiving_country,
    corridor      = paste(origin, "\u2192", destination),
    corridor_name = corridor_name,
    Model         = "Observed",
    mean_value    = pred_q50,
    lower         = pred_q05,
    upper         = pred_q95
  )

# Simulated corridor summaries (Gravity & MST)
sim_corr <- chosen_corridors %>% 
  rename(mean_value = mean_sim,
         lower = min_sim,
         upper = max_sim,
         Model = model) %>% 
  mutate(Year = as.numeric(year))


# Combine observed & simulated data for plotting 
ribbons_df <- bind_rows(
  sim_corr,
  obs_corr
) %>%
  mutate(
    Model = factor(Model, levels = c("Observed", "Gravity", "MST"))
  )

lines_df <- ribbons_df %>%
  transmute(
    Year,
    corridor_name,
    Model,
    flow = mean_value
  ) %>% 
  mutate(
    Model = factor(Model, levels = c("Observed", "Gravity", "MST"))
  )

# Aesthetics
corr_colors <- c(
  "Observed" = "#000000",
  "Gravity"  = "#0072B2",
  "MST"      = "#E69F00"
)

corr_fills <- c(
  "Observed" = "#3B3B3B", 
  "Gravity"  = "#A6CEE3",
  "MST"      = "#FDBF6F"
)

#corr_colors <- c(
#  "Observed" = "#000000",
#  "Gravity"  = "#009E73" ,
#  "MST"      = "#CC79A7"
#)

#corr_fills <- c(
#  "Observed" = "#3B3B3B", 
#  "Gravity"  = "#66D1AB",
#  "MST"      = "#E0AFCF"
#)


corr_linetypes <- c(
  "Observed" = "solid",
  "Gravity"  = "dotted",
  "MST"      = "longdash"
)


p_corridors <- ggplot() +
  facet_wrap(
    ~ corridor_name,
    scales = "free_y",
    ncol   = 2,
    strip.position = "top"
  ) +
  (
    geom_ribbon(
      data = ribbons_df,
      aes(
        x     = Year,
        ymin  = lower / 1000,
        ymax  = upper / 1000,
        fill  = Model,
        colour = Model,
        linetype = Model
      ),
      alpha     = 0.35,
      linewidth = 0.5
    ) %>%
      partition(vars(Model)) %>%   # blend per model
      blend("multiply")
  ) +
  geom_line(
    data = lines_df,
    aes(
      x        = Year,
      y        = flow / 1000,
      color    = Model,
      linetype = Model
    ),
    linewidth = 1.0
  ) +
  scale_color_manual(values   = corr_colors) +
  scale_fill_manual(values    = corr_fills) +
  scale_linetype_manual(values = corr_linetypes) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4),
    labels = scales::label_number(accuracy = 1)
  ) +
  theme_ready_pub(base_size = 18) +
  labs(
    title = "Observed and Simulated Flows in Selected Corridors, 2002–2021",
    subtitle = paste0(
      "Observed values in black, with median values (lines) and interquartile ranges (ribbons).\n",
      "Simulated values in blue (gravity) and orange (MST), with simulation means (lines)\n",
      "and entire range of distributions (ribbons). Y-axis scales differ by country."
    ),
    x        = "Year",
    y        = "Migration volume (thousands)",
    fill     = "Model:",
    color    = "Model:",
    linetype = "Model:"
  ) +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.key.width  = unit(2.5, "cm"),
    legend.key.height = unit(0.7, "cm")
  )

p_corridors 

ggsave(
  p_corridors,
  file   = here("figures", "selected_corridors_flows.png"),
  device = png,
  type   = "cairo",
  width  = 13,
  height = 13,
  dpi    = 300
)

# 6. Scatterplots of simulated vs observed by region ---------------------------


# Shapes by model
model_shapes <- c(
  "Gravity" = 1,  # hollow circle
  "MST"     = 3   # cross
)

region_map <- c(
  "Eastern Europe"  = "EE",
  "Northern Europe" = "NE",
  "Western Europe"  = "WE",
  "Southern Europe" = "SE"
)

corridor_fit <- corridor_fit %>%
  mutate(
    region_corr = region_pair %>%
      str_replace_all(region_map) %>%
      str_replace(" → ", " → ")
  )

p_sim_vs_obs_region_gravity <- ggplot(
  corridor_fit %>% 
    filter(model == "Gravity"),
  aes(
    x     = pred_q50,
    y     = mean_sim,
    shape = model,
    colour = year
  )
) +
  facet_wrap(~ region_corr, scales = "free") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_point(
    alpha = 0.4,
    size  = 2
  ) +
  scale_shape_manual(values = model_shapes) +
  scale_color_viridis_c(option = "viridis", end = 0.95) +
  scale_x_continuous(
    trans  = "log10",
    labels = label_number(accuracy = 1),
    breaks = 10^(1:4)
  ) +
  scale_y_continuous(
    trans  = "log10",
    labels = label_number(accuracy = 1),
    breaks = 10^(1:4)
  ) +
  theme_ready_pub(base_size = 18) +
  theme(axis.text = element_text(size = rel(0.70))) +
  labs(
    x      = "Observed flows (log scale)",
    y      = "Simulated flows (mean across simulations, log scale)",
    colour = "Year:",
    shape  = "Model:",
    title  = "Gravity: Simulated vs Observed Flows by Corridor–Year and Region Pair",
    subtitle = "Each point is a corridor–year observation; dashed lines indicate perfect simulation.",
    caption = "Region labels: EE = Eastern Europe, NE = Northern Europe, SE = Southern Europe, WE = Western Europe."
  ) +
  guides(
    colour = guide_colorbar(barwidth = unit(8, "cm"),
                            barheight = unit(1.5, "cm")),
    shape  = guide_legend(override.aes = list(size = 6, alpha = 1))
  )

ggsave(
  p_sim_vs_obs_region_gravity,
  file   = here("figures", "scatter_sim_vs_obs_by_region_gravity.png"),
  device = png,
  type   = "cairo",
  width  = 14.5,
  height = 15,
  units  = "in",
  dpi    = 300
)
 


p_sim_vs_obs_region_mst <- ggplot(
  corridor_fit %>% 
    filter(model == "MST"),
  aes(
    x     = pred_q50,
    y     = mean_sim,
    shape = model,
    colour = year
  )
) +
  facet_wrap(~ region_corr, scales = "free") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_point(
    alpha = 0.4,
    size  = 2
  ) +
  scale_shape_manual(values = model_shapes) +
  scale_color_viridis_c(option = "viridis", end = 0.95) +
  scale_x_continuous(
    trans  = "log10",
    labels = label_number(accuracy = 1),
    breaks = 10^(1:4)
  ) +
  scale_y_continuous(
    trans  = "log10",
    labels = label_number(accuracy = 1),
    breaks = 10^(1:4)
  ) +
  theme_ready_pub(base_size = 18) +
  theme(axis.text = element_text(size = rel(0.70))) +
  labs(
    x      = "Observed flows (log scale)",
    y      = "Simulated flows (mean across simulations, log scale)",
    colour = "Year:",
    shape  = "Model:",
    title  = "MST: Simulated vs Observed Flows by Corridor–Year and Region Pair",
    subtitle = "Each point is a corridor–year observation; dashed lines indicate perfect simulation.",
    caption = "Region labels: EE = Eastern Europe, NE = Northern Europe, SE = Southern Europe, WE = Western Europe."
  ) +
  guides(
    colour = guide_colorbar(barwidth = unit(8, "cm"),
                            barheight = unit(1.5, "cm")),
    shape  = guide_legend(override.aes = list(size = 6, alpha = 1))
  )

ggsave(
  p_sim_vs_obs_region_mst,
  file   = here("figures", "scatter_sim_vs_obs_by_region_mst.png"),
  device = png,
  type   = "cairo",
  width  = 14.5,
  height = 15,
  units  = "in",
  dpi    = 300
)
