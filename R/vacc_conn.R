### Functions for running most pared down IBM possible
## Malavika Rajeev
## December 2018
rm(list = ls())
library(doMPI)
cl <- startMPIcluster()
clusterSize(cl) # this just tells you how many you've got
registerDoMPI(cl)

## Libraries
library(tidyverse)
library(magrittr)
library(lubridate)
library(rgdal)
library(raster)
library(ISOweek)
library(data.table)
library(foreach)

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
y1 <- 2002
ylast <- 2015
steps <- 52
tmax <- (ylast - y1) * steps
start_date <- "01-01-2002"

## Get vacc campaigns at vill level
campaigns <- get.campaigns.WM(vacc = vacc_data, pop = pop_data, shape = SD_shape, threshold = 7)
campaigns$week <- get.consec(campaigns$date_med, format_date = "%Y-%m-%d", start = "01-01-2002", 
                             format_start = "%d-%m-%Y", year1 = 2002, tstep = "week", get.info = FALSE)
# 
# campaigns %>%
#   group_by(villcode) %>%
#   arrange(week) %>%
#   mutate(weeks_between = week - lag(week)) -> check

ts <- as.data.frame(1:tmax) 
names(ts) <- "week"
ts %>% 
  left_join(campaigns) %>%
  dplyr::select(week, total, villcode) %>%
  spread(week, total, fill = 0) %>%
  filter(!is.na(villcode)) -> vacc_mat
vacc_mat <- as.data.table(list(villcode = grid_data$villcode, 
                               cell_id = grid_data$cell_id))[vacc_mat, on = "villcode"]
vacc_mat <- vacc_mat[order(match(vacc_mat$cell_id, grid_data$cell_id)), ] 
vacc_mat <- as.matrix(vacc_mat[, c("cell_id","villcode"):=NULL]) ## getting rid of id cols

## rep incursions (so as to include scaling factor)
iota = 1; scale_iota = 1;
incs <- ifelse(length(scale_iota) > 1, iota*scale_iota, rep(iota, length(tmax)))

##' Sim cov ----------------------------------------------------------------------------------------
##' ## Need for time step functions within data cleaning
burn_in <- 2
sim_yrs <- 8
tmax <- (burn_in + sim_yrs) * 52
burn_past <- tmax - burn_in*52

## Get vacc campaigns at vill level
## rep incursions (so as to include scaling factor)
iota = 1; scale_iota = 1;
incs <- ifelse(length(scale_iota) > 1, iota*scale_iota, rep(iota, length(tmax)))
perc_val <- seq(0, 1, by = 0.1)
cov_val <- seq(0, 1, by = 0.1)
nsim <- 1000

system.time({
  perc_df <- foreach(i = 1:length(perc_val), .combine = rbind) %:%
    foreach(j = 1:length(cov_val), .combine = rbind) %:%
    foreach(k = 1:nsim, .combine = rbind, 
            .packages = c("igraph", "data.table", "dplyr", "sp", "raster")) %dopar% {
              
              vacc_mat <- sim.campaigns(data = grid_data, vills = unique(grid_data$villcode),
                                        sim_years = sim_yrs, burn_in_years = burn_in, 
                                        vill_weeks = sample(1:75, 75, replace = TRUE))
              
              check <- sim.IBM(R_0 = 0, k = 0.4, inc_rate = 1, p_revacc = 1, I_seed = 10,
                               p_vacc = 1, start_vacc = 0, perc_val = perc_val[i],
                               vacc_sim = "sim_vill", cov_val = cov_val[j],
                               return_coords = FALSE)
              ## Get cov + conn metric across time series (covered vs. not) (+ @ 0.2 & 0.3)

            }
})

write.csv(perc_conn, "output/perc_conn.csv")

### Then just close it out at the end
print("Done:)")
closeCluster(cl) ## for doMPI
mpi.quit()
