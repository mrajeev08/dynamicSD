## Gif potentially
get.gif <- function (coords = I_coords, shape = SD_shape, tstep = tmax, gif_name,
                     nweeks = 4){
  saveGIF ({ani.options(convert = '/opt/local/bin/convert', interval = 0.5, nmax = 50)
    for(t in seq(1, tmax, by = nweeks)){
      cases <- coords[tstep <= t & tstep >= t - 1]
      plot(shape)
      points(jitter(cases$x_coord, factor = 5), jitter(cases$y_coord, factor = 5), 
             col = ifelse(cases$progen_ID == 0, "pink", "red"), pch = 1, cex = 2)
      progens <- coords[progen_ID %in% cases$ID]
      points(progens$x_coord, progens$y_coord, col = alpha("blue", 0.2), pch = 20, 
             cex = 2)
      ani.pause()
    }
  },
  movie.name = gif_name, ani.width = 600, ani.height = 600, clean=TRUE)
}
get.gif(gif_name = "test.gif")

get.gif.vacc <- function(rast_SD = SD_raster, cov = cov_mat, tstep = tmax, 
                         nweeks = 4, gif_name){
  saveGIF ({ani.options(convert = '/opt/local/bin/convert', interval = 0.5, nmax = 50)
    for(t in seq(1, tmax, by = nweeks)){
      rast <- rast_SD
      values(rast) <-  1 - S[row_id[match(values(rast), cell_id)], t]/
        N[row_id[match(values(rast), cell_id)], t]  
      brks <- seq(0, 1, by = 0.1)
      pal <- colorRampPalette(c("white", "blue"))(length(brks))
      plot(rast, axes = FALSE, box = FALSE, breaks = brks, col = pal)
      ani.pause()
    }
  },
  movie.name = gif_name, ani.width = 600, ani.height = 600, clean=TRUE)
}
get.gif.vacc(gif_name = "testvacc.gif")

get.gif.combined <- function(coords = I_coords, 
                             rast_SD = SD_raster, cov = cov_mat, tstep = tmax,
                             nweeks = 4, gif_name){
  saveGIF ({ani.options(convert = '/opt/local/bin/convert', interval = 0.5, nmax = 50)
    for(t in seq(1, tmax, by = nweeks)){
      rast <- rast_SD
      values(rast) <-  1 - S[row_id[match(values(rast), cell_id)], t]/
        N[row_id[match(values(rast), cell_id)], t]  
      brks <- seq(0, 1, by = 0.1)
      pal <- colorRampPalette(c("white", "blue"))(length(brks))
      plot(rast, axes = FALSE, box = FALSE, breaks = brks, col = pal, colNA = "grey")
  
      ## cases
      cases <- coords[tstep <= t & tstep >= t - 1]
      points(jitter(cases$x_coord, factor = 5), jitter(cases$y_coord, factor = 5), 
             col = ifelse(cases$progen_ID == 0, "darkred", "red"), pch = 1, cex = 2)
      progens <- coords[progen_ID %in% cases$ID]
      points(progens$x_coord, progens$y_coord, col = alpha("orange", 0.2), pch = 20, 
             cex = 2)
      ani.pause()
    }
  },
  movie.name = gif_name, ani.width = 600, ani.height = 600, clean=TRUE)
}
get.gif.combined(gif_name = "testboth.gif")
