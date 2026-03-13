# scripts/utils/plotting_helpers.R

library(ggthemes)
library(ggplot2)
library(grid)
library(lemon)


# Legend shifters
shift_legend <- function(p) {
  if (!inherits(p, "gtable")) {
    if (inherits(p, "ggplot")) {
      gp <- ggplotGrob(p)
    } else {
      message("Not a ggplot object; returning original.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(
    facet.panels,
    function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]),
    USE.NAMES = FALSE
  )
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  if (length(empty.facet.panels) == 0) {
    message("No empty facet panels to place legend into; returning original.")
    return(p)
  }
  
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  
  reposition_legend(p, "center", panel = names, plot = FALSE)
}

shift_legend2 <- function(p) {
  if (!inherits(p, "gtable")) {
    if (inherits(p, "ggplot")) {
      gp <- ggplotGrob(p)
    } else {
      message("Not a ggplot object; returning original.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(
    facet.panels,
    function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]),
    USE.NAMES = FALSE
  )
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  if (length(empty.facet.panels) == 0) {
    message("No empty facet panels to place legend into; returning original.")
    return(p)
  }
  
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  
  reposition_legend(p, "center", panel = names)
}

# Publication-ready theme
# Publication-ready theme
theme_ready_pub <- function(base_size = 18, base_family = "Public Sans") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face   = "bold",
        size   = rel(1.5),
        color  = "#1A242F",
        hjust  = 0.5,
        margin = margin(t = 0, r = 0, b = 15, l = 0)
      ),
      plot.subtitle = element_text(
        hjust  = 0.5,
        size   = rel(1.2),
        color  = "#1A242F",
        margin = margin(t = 0, r = 0, b = 25, l = 0)
      ),
      plot.caption = element_text(
        hjust  = 1,
        color  = "#474F58",
        margin = margin(t = 15, r = 5, b = 0, l = -5)
      ),
      text             = element_text(color = "#1A242F"),
      panel.background = element_rect(colour = NA, fill = "white"),
      plot.background  = element_rect(colour = NA, fill = "white"),
      panel.border     = element_rect(colour = NA),
      panel.spacing    = unit(2, "lines"),
      axis.title       = element_text(face = "bold", size = rel(1.1)),
      axis.title.y     = element_text(
        angle  = 90,
        colour = "#1A242F",
        vjust  = 2,
        hjust  = 1
      ),
      axis.title.x     = element_text(
        vjust  = -0.2,
        colour = "#1A242F",
        hjust  = 0
      ),
      axis.text        = element_text(colour = "#474F58", size = rel(0.9)),
      axis.line        = element_line(colour = "#474F58"),
      axis.ticks       = element_line(colour = "#474F58"),
      panel.grid.major = element_line(colour = "#E3E4E5", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      
      # Legend
      legend.key       = element_rect(colour = NA, fill = NA),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.key.size  = unit(0.6, "cm"),
      legend.title     = element_text(face = "bold", color = "#1A242F"),
      legend.text      = element_text(size = rel(1)),
      plot.margin      = unit(c(10, 20, 10, 20), "pt"),
      
      # Facet strip with NO BOX — only text
      strip.background = element_blank(),
      strip.text       = element_text(
        face = "bold",
        color = "#1A242F",
        size = rel(1.0)
      )
    )
}