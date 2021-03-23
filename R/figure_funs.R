# figure funs: theming

plotpars <- function() {
  
  theme_proj <- function() {
    cowplot::theme_minimal_hgrid(font_size = 12) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  theme_map <- function() {
    cowplot::theme_map(font_size = 12)
  }
  
  
  summ_stat_labs <- read_csv("analysis/out/summ_stat_labs.csv")
  
  mod_labs <- read_csv("analysis/out/model_labs.csv")
  
  type_cols <-
    tribble(~color, ~type,
            "#009E73", "Spatiotemporal",
            "#D55E00", "Spatial",
            "#0072B2", "Temporal",
            "#000000", "Model based")
  
  summ_stat_labs %<>%
    left_join(type_cols) %>%
    mutate(lab = coalesce(lab, summstat),
           name = glue("<p style='color:{color}'>{lab}</p>"))
  
  # summ stat labs
  slabs <- summ_stat_labs$name
  names(slabs) <- summ_stat_labs$summstat
  
  # summ stat cols
  scols <- type_cols$color
  names(scols) <- type_cols$type
  scale_labs <- c("1x1 km", "Village")
  names(scale_labs) <- c("grid", "vill")
  
  # model run labs (i.e. full or subsampled)
  run_labs <- c("Full simulation set", "Subsampled set (75%, N = 3)")
  names(run_labs) <- c("full", "se")
  fct_incs <- c("Introductions \n estimated", "Introductions \n hardwired")
  names(fct_incs) <- c("Estimated", "Fixed")
  
  as.list.environment(environment())
}

# somewhat hackish solution to:
# https://twitter.com/EamonCaddigan/status/646759751242620928
# based mostly on copy/pasting from ggplot2 geom_violin source:
# https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r

library(ggplot2)
library(dplyr)


"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )


