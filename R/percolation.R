### Functions for running most pared down IBM possible
## Malavika Rajeev
## December 2018
rm(list = ls())

## Libraries
library(tidyverse)
library(magrittr)
library(lubridate)
library(rgdal)
library(raster)
library(ISOweek)
library(data.table)
library(foreach)
library(igraph)

## Source in functions
source("R/data_functions.R")
source("R/utils.R")
source("R/sim_functions.R")
source("R/sim.IBM.R")
source("R/vacc_functions.R")

## Data
vacc_data <- read_csv(paste0("data/", list.files("data/", pattern = "Vaccination")))
SD_shape <- readOGR("data/SD_shape/Serengeti_villages UTM_region.shp")
census_data <- read.csv("data/SDcompiled.csv")
pop_data <- read.csv("data/SerengetiPop.csv")

## Get raster data
r <- raster(SD_shape)
res(r) <- 1000
SD_raster <- rasterize(SD_shape, r)

## Data by cellID
grid_data <- get.demdata(shapefile = SD_shape, res_meters = 1000, pop = pop_data, census = census_data)
values(SD_raster) <- grid_data$cell_id ## cell IDs as values

## create coordinates of raster
coords <- coordinates(SD_raster)
grid_data$x_coord <- coords[, 1]
grid_data$y_coord <- coords[, 2]

## Filter to ones inside SD villages 
grid_data <- as.data.table(subset(grid_data, !is.na(village_ID) & !is.na(start_pop)))
# subset out ones not in SD vills and also change to data table to speed things up!
## also subsetting out empty cells...need to think if I NEED to populate these somehow...

## Need for time step functions within data cleaning
tmax <- 15 * 52

## Get vacc campaigns at vill level
## rep incursions (so as to include scaling factor)
iota = 1; scale_iota = 1;
incs <- ifelse(length(scale_iota) > 1, iota*scale_iota, rep(iota, length(tmax)))
perc_val <- seq(0, 1, by = 0.5)
cov_val <- seq(0, 1, by = 0.5)
nsim = 2
burn_past <- tmax - 5*52

sum_times <- function(vector, steps, na.rm=TRUE) {    # 'matrix'
  nv <- length(vector)
  if (nv %% steps)
    vector[ceiling(nv / steps) * steps] <- NA
  colSums(matrix(vector, steps), na.rm = na.rm)
}

multi_cbind <- function(x, ...) {  
  mapply(cbind, x, ..., SIMPLIFY = FALSE)
}

multi_rbind <- function(x, ...) {  
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}  

perc_df <- foreach(i = 1:length(perc_val), .combine = rbind) %:%
  foreach(j = 1:length(cov_val), .combine = rbind) %:%
    foreach(k = 1:nsim, .combine = rbind) %dopar% {
      
      vacc_mat <- sim.campaigns(data = grid_data, vills = unique(grid_data$villcode),
                                sim_years = 10, burn_in_years = 5, 
                                vill_weeks = sample(1:75, 75, replace = TRUE))
      
      check <- sim.IBM(R_0 = 1.1, k = 0.4, inc_rate = 1, p_revacc = 1, I_seed = 10,
                     p_vacc = 1, start_vacc = 0, perc_val = 1,
                     vacc_sim = "sim_cell", cov_val = 1,
                     return_coords = TRUE)
      
      ## Get number of cases seeded by a local case in last half of sim period
      I_coords <- check[["I_coords"]]
      I_coords$second_progen <- I_coords$progen_ID[match(I_coords$progen_ID, I_coords$ID)]
      secondaries <- nrow(I_coords %>% filter(second_progen != 0, !is.na(second_progen), tstep >= burn_past))
      secondaries <- as.data.frame(list(perc = perc_val[i], cov = cov_val[j], sim = k, 
                                        secondaries = secondaries))
      
      I_coords %>%
        dplyr::select(from = progen_ID, to = ID) %>%
        filter(from != 0) -> edge_df
      gcheck <- graph_from_data_frame(edge_df, directed = TRUE)
      hist(components(gcheck)$csize)
      chains <- as.data.frame(list(case_ID = as.numeric(names(components(gcheck)$membership)), 
                                   chain_ID = components(gcheck)$membership))
      I_coords %>%
        left_join(chains, by = c("ID" = "case_ID")) %>%
        mutate(chain_ID = if_else(is.na(chain_ID), 0, chain_ID), 
               order_id = 1:n()) %>%
        group_by(chain_ID) %>%
        filter(tstep >= burn_past, chain_ID != 0) %>%
        summarize(min = min(tstep), max = max(tstep), 
                  length = max - min) -> chains
      
      as.data.frame(list(perc = perc_val[i], cov = cov_val[j], sim = k,
                         chain_length = chains$length, secondaries = secondaries))
    }

write.csv(perc_df, "output/perc_df.csv")