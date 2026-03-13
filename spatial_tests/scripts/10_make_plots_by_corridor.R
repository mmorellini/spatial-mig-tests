library(dplyr)
library(ggplot2)
library(ggblend)
library(here)



corridors_lab <- sim_summaries %>%
  arrange(desc(mean_sim)) %>%      # highest overall_mean first
  mutate(
    corridor = paste(origin, "\u2192", destination)
  ) %>%
  mutate(corridor = factor(corridor))


obs_corr <- merged %>%
  semi_join(
    sim_summaries %>% select(origin, destination),
    by = c("sending_country" = "origin",
           "receiving_country" = "destination")
  ) %>%
  transmute(
    Year      = as.numeric(year),
    origin    = sending_country,
    destination = receiving_country,
    corridor  = paste(origin, "\u2192", destination),
    Model     = "Observed",
    mean_value = pred_q50,
    lower      = pred_q05,
    upper      = pred_q95
  ) %>%
  # keep only our 31 corridors and impose facet order
  semi_join(corridors_lab, by = c("origin", "destination"))

## 3. Simulated corridor summaries (Gravity & MST) -----------------------------

sim_corr <- sim_summaries %>% 
  rename(mean_value = mean_sim,
         lower = min_sim,
         upper = max_sim,
         Model = model) %>% 
  mutate(Year = as.numeric(year),
         corridor = paste(origin, "\u2192", destination)
  ) %>%
  mutate(corridor = factor(corridor))


## 4. Combine observed & simulated data for plotting ---------------------------

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
    corridor,
    Model,
    flow = mean_value
  ) %>% 
  mutate(
    Model = factor(Model, levels = c("Observed", "Gravity", "MST"))
  )

corridor_dir <- here("figures", "corridors")
if (!dir.exists(corridor_dir)) {
  dir.create(corridor_dir, recursive = TRUE)
}

corridors_all <- sort(unique(ribbons_df$corridor))


## Aesthetics for corridors

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

corr_linetypes <- c(
  "Observed" = "solid",
  "Gravity"  = "dotted",
  "MST"      = "longdash"
)


for (corr in corridors_all) {
  message("Plotting corridor: ", corr)
  
  # Filter data for this corridor only
  ribbons_corr <- ribbons_df %>% filter(corridor == corr)
  lines_corr   <- lines_df   %>% filter(corridor == corr)
  
  # Make a safe filename, e.g. "AT_DE" instead of "AT → DE"
  fname_safe <- gsub("[^A-Za-z0-9]+", "_", corr)
  
  p_corr <- ggplot() +
    # BLENDED RIBBONS for this corridor
    (
      geom_ribbon(
        data = ribbons_corr,
        aes(
          x     = Year,
          ymin  = lower,   # keep your scaling
          ymax  = upper ,
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
    # LINES (mean values) for this corridor
    geom_line(
      data = lines_corr,
      aes(
        x        = Year,
        y        = flow,
        color    = Model,
        linetype = Model
      ),
      linewidth = 1.0
    ) +
    scale_color_manual(values   = corr_colors) +
    scale_fill_manual(values    = corr_fills) +
    scale_linetype_manual(values = corr_linetypes) +
    scale_x_continuous(breaks = seq(2002, 2021, 3)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 4),
      labels = scales::label_number(accuracy = 1)
    ) +
    theme_ready_pub(base_size = 18) +
    labs(
      title = paste0("Observed and Simulated Flows, ", corr, ", 2002–2021"),
      subtitle = paste0(
        "Observed Bayesian flows (black) versus simulations from Gravity (blue) and MST (yellow).\n",
        "Lines report the mean of each distribution; ribbons show min–max across simulations (5–95% for observed).\n",
        "Y-axis in hundreds."
      ),
      x        = "Year",
      y        = "Migration volume (hundreds)",
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
  
  
  ggsave(
    filename = file.path(corridor_dir, paste0("corridor_", fname_safe, ".png")),
    plot     = p_corr,
    device   = png,
    type     = "cairo",
    width    = 16,
    height   = 12,
    dpi      = 300
  )
}
