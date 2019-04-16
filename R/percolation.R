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
tmax <- 20 * 52

## Get vacc campaigns at vill level
vacc_mat <- sim.campaigns(data = grid_data, vills = unique(grid_data$villcode),
                          sim_years = 15, burn_in_years = 5, 
                          vill_weeks = sample(1:75, 75, replace = TRUE))
  
## rep incursions (so as to include scaling factor)
iota = 1; scale_iota = 1;
incs <- ifelse(length(scale_iota) > 1, iota*scale_iota, rep(iota, length(tmax)))

check <- sim.IBM(R_0 = 1.1, k = 0.5, inc_rate = 1, p_revacc = 1, 
                 p_vacc = 1, start_vacc = 0, perc_val = 0.50,
                 vacc_sim = "sim_cell", cov_val = 1,
                 return_coords = TRUE)

## Get number of cases seeded by a local case in last half of sim period
I_coords <- check[["I_coords"]]
I_coords$second_progen <- I_coords$progen_ID[match(I_coords$progen_ID, I_coords$ID)]
nrow(I_coords %>% filter(second_progen != 0, !is.na(second_progen), tstep >= tmax - 15*52/2))

sum_times <- function(vector, steps, na.rm=TRUE) {    # 'matrix'
  nv <- length(vector)
  if (nv %% steps)
    vector[ceiling(nv / steps) * steps] <- NA
  colSums(matrix(vector, steps), na.rm = na.rm)
}

N <- check[["N"]]
S <- check[["S"]]
E <- check[["E"]]
I_all <- check[["I_all"]]
I_dist <- check[["I_dist"]]

mIobs <- apply(I_dist[,1:(ncol(I_dist)-3)], 1 , sum_times, steps = 4)
plot(rowSums(mIobs, na.rm = TRUE), type = "l")
plot(colSums(I_all), col = "red", type = "l")
plot(colSums(S)/colSums(N), col = "blue", type = "l", ylim = c(0, 1))
plot(1 - colSums(S)/colSums(N), col = "blue", type = "l", ylim = c(0, 1)) 
max(1 - colSums(S)/colSums(N))
mean(1 - colSums(S)/colSums(N))

plot(colSums(N), col = "blue", type = "l")

mNobs <- N[, seq(1, tmax, by = 4)]

