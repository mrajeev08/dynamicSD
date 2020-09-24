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
burn_in <- 5
sim_yrs <- 10
tmax <- (burn_in + sim_yrs) * 52
burn_past <- tmax - burn_in*52

## Get vacc campaigns at vill level
## rep incursions (so as to include scaling factor)
perc_val <- seq(0, 1, by = 0.1)
cov_val <- seq(0, 1, by = 0.1)
nsim <- 500

## Make directory if not already one
dir_name <- paste0("output/perc")

if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

##' Village level sims with heterogeneity from the census
##' ------------------------------------------------------------------------------------------------
system.time({
  perc_df <- foreach(i = 1:length(perc_val), .combine = rbind) %:%
  foreach(j = 1:length(cov_val), .combine = rbind) %:%
    foreach(k = 1:nsim, .combine = rbind, 
            .packages = c("igraph", "data.table", "dplyr", "sp", "raster")) %dopar% {
      
      vacc_mat <- sim.campaigns(data = grid_data, vills = unique(grid_data$villcode),
                                sim_years = sim_yrs, burn_in_years = burn_in, 
                                vill_weeks = sample(1:52, 75, replace = TRUE))
      
      check <- sim.IBM(R_0 = 1.1, k = 0.5, inc_rate = 1, p_revacc = 1, I_seed = 1,
                     p_vacc = grid_data$cov, start_vacc = 0, perc_val = perc_val[i],
                     vacc_sim = "sim_vill", cov_val = cov_val[j],
                     return_coords = TRUE)
      
      ## Get number of cases seeded by a local case in last half of sim period
      I_coords <- check[["I_coords"]]
      I_coords$second_progen <- I_coords$progen_ID[match(I_coords$progen_ID, I_coords$ID)]
      secondaries <- nrow(I_coords %>% filter(second_progen != 0, !is.na(second_progen), tstep >= burn_past))
      
      I_coords %>%
        dplyr::select(from = progen_ID, to = ID) %>%
        filter(from != 0) -> edge_df
      gcheck <- graph_from_data_frame(edge_df, directed = TRUE)
      chains <- as.data.frame(list(case_ID = as.numeric(names(components(gcheck)$membership)), 
                                   chain_ID = components(gcheck)$membership))
      I_coords %>%
        left_join(chains, by = c("ID" = "case_ID")) %>%
        mutate(chain_ID = if_else(is.na(chain_ID), 0, chain_ID), 
               order_id = 1:n(), 
               step_id = paste0(chain_ID, "_", progen_ID)) %>%
        group_by(chain_ID) %>%
        mutate(step_id = as.numeric(as.factor(step_id))) %>%
        filter(tstep >= burn_past, chain_ID != 0) %>% ## only those past the burn in and vacc periods
        summarize(chain_size = max(step_id)) -> chains
      
      perc_df <- as.data.frame(list(perc = perc_val[i], cov = cov_val[j], sim = k,
                         chain_length = max(chains$chain_size), secondaries = secondaries))
    }
})
write.csv(perc_df, "output/perc/perc_vill.csv")


##' Village level sims with scaled down incursions as disrict cov decreases
##' ------------------------------------------------------------------------------------------------
system.time({
  perc_df <- foreach(i = 1:length(perc_val), .combine = rbind) %:%
    foreach(j = 1:length(cov_val), .combine = rbind) %:%
    foreach(k = 1:nsim, .combine = rbind, 
            .packages = c("igraph", "data.table", "dplyr", "sp", "raster")) %dopar% {
              
              vacc_mat <- sim.campaigns(data = grid_data, vills = unique(grid_data$villcode),
                                        sim_years = sim_yrs, burn_in_years = burn_in, 
                                        vill_weeks = sample(1:52, 75, replace = TRUE))
              
              check <- sim.IBM(R_0 = 1.1, k = 0.5, inc_rate = 1, p_revacc = 1, I_seed = 1,
                               p_vacc = grid_data$cov, start_vacc = 0, perc_val = perc_val[i],
                               vacc_sim = "sim_vill", cov_val = cov_val[j],
                               return_coords = TRUE, scale_incs = TRUE)
              
              ## Get number of cases seeded by a local case in last half of sim period
              I_coords <- check[["I_coords"]]
              I_coords$second_progen <- I_coords$progen_ID[match(I_coords$progen_ID, I_coords$ID)]
              secondaries <- nrow(I_coords %>% filter(second_progen != 0, !is.na(second_progen), tstep >= burn_past))
              
              I_coords %>%
                dplyr::select(from = progen_ID, to = ID) %>%
                filter(from != 0) -> edge_df
              gcheck <- graph_from_data_frame(edge_df, directed = TRUE)
              chains <- as.data.frame(list(case_ID = as.numeric(names(components(gcheck)$membership)), 
                                           chain_ID = components(gcheck)$membership))
              I_coords %>%
                left_join(chains, by = c("ID" = "case_ID")) %>%
                mutate(chain_ID = if_else(is.na(chain_ID), 0, chain_ID), 
                       order_id = 1:n(), 
                       step_id = paste0(chain_ID, "_", progen_ID)) %>%
                group_by(chain_ID) %>%
                mutate(step_id = as.numeric(as.factor(step_id))) %>%
                filter(tstep >= burn_past, chain_ID != 0) %>% ## only those past the burn in and vacc periods
                summarize(chain_size = max(step_id)) -> chains
              
              perc_df <- as.data.frame(list(perc = perc_val[i], cov = cov_val[j], sim = k,
                                            chain_length = max(chains$chain_size), 
                                            secondaries = secondaries))
            }
})
write.csv(perc_df, "output/perc/perc_vill_scaleIntros.csv")

##' Village level sims with pups vaccinated at 25%
##' ------------------------------------------------------------------------------------------------
pup_vacc <- c(0.25, 0.5, 0.75)

system.time({
  perc_df <- foreach(i = 1:length(perc_val), .combine = rbind) %:%
    foreach(j = 1:length(cov_val), .combine = rbind) %:%
    foreach(k = 1:length(nsim), .combine = rbind) %:%
    foreach(l = 1:length(pup_vacc), .combine = rbind, 
            .packages = c("igraph", "data.table", "dplyr", "sp", "raster")) %dopar% {
              
              vacc_mat <- sim.campaigns(data = grid_data, vills = unique(grid_data$villcode),
                                        sim_years = sim_yrs, burn_in_years = burn_in, 
                                        vill_weeks = sample(1:52, 75, replace = TRUE))
              
              check <- sim.IBM(R_0 = 1.1, k = 0.5, inc_rate = 1, p_revacc = 1, I_seed = 1,
                               p_vacc = grid_data$cov, start_vacc = 0, perc_val = perc_val[i],
                               vacc_sim = "sim_vill", cov_val = cov_val[j],
                               return_coords = TRUE, 
                               pup_vacc = c(rep(0, burn_in*52), rep(pup_vacc[l], sim_yrs*52)))
              
              ## Get number of cases seeded by a local case in last half of sim period
              I_coords <- check[["I_coords"]]
              I_coords$second_progen <- I_coords$progen_ID[match(I_coords$progen_ID, I_coords$ID)]
              secondaries <- nrow(I_coords %>% filter(second_progen != 0, !is.na(second_progen), 
                                                      tstep >= burn_past))
              
              I_coords %>%
                dplyr::select(from = progen_ID, to = ID) %>%
                filter(from != 0) -> edge_df
              gcheck <- graph_from_data_frame(edge_df, directed = TRUE)
              chains <- as.data.frame(list(case_ID = as.numeric(names(components(gcheck)$membership)), 
                                           chain_ID = components(gcheck)$membership))
              I_coords %>%
                left_join(chains, by = c("ID" = "case_ID")) %>%
                mutate(chain_ID = if_else(is.na(chain_ID), 0, chain_ID), 
                       order_id = 1:n(), 
                       step_id = paste0(chain_ID, "_", progen_ID)) %>%
                group_by(chain_ID) %>%
                mutate(step_id = as.numeric(as.factor(step_id))) %>%
                filter(tstep >= burn_past, chain_ID != 0) %>% ## only those past the burn in and vacc periods
                summarize(chain_size = max(step_id)) -> chains
              
              perc_df <- as.data.frame(list(perc = perc_val[i], cov = cov_val[j], sim = k,
                                            chain_length = max(chains$chain_size), 
                                            secondaries = secondaries, pup_vacc = pup_vacc[l]))
            }
})
write.csv(perc_df, "output/perc/perc_vill_vaccPups.csv")

##' Flat vaccination implemented to pup cov
##'  ------------------------------------------------------------------------------------------------
perc_val <- seq(0, 1, by = 0.1)
pup_vacc <- seq(0, 1, by = 0.1)
vacc_mat <- matrix
system.time({
  perc_df <- foreach(i = 1:length(perc_val), .combine = rbind) %:%
    foreach(k = 1:length(nsim), .combine = rbind) %:%
    foreach(l = 1:length(pup_vacc), .combine = rbind, 
            .packages = c("igraph", "data.table", "dplyr", "sp", "raster")) %dopar% {
              vacc_mat <- matrix(0, nrow = nrow(grid_data), ncol = tmax)
              check <- sim.IBM(R_0 = 1.1, k = 0.5, inc_rate = 1, p_revacc = 1, I_seed = 1,
                               p_vacc = grid_data$cov, start_vacc = 0, perc_val = perc_val[i],
                               vacc_sim = "sim_vill", cov_val = cov_val[j],
                               return_coords = TRUE, 
                               pup_vacc = c(rep(0, burn_in*52), rep(pup_vacc[l], sim_yrs*52)))
              
              ## Get number of cases seeded by a local case in last half of sim period
              I_coords <- check[["I_coords"]]
              I_coords$second_progen <- I_coords$progen_ID[match(I_coords$progen_ID, I_coords$ID)]
              secondaries <- nrow(I_coords %>% filter(second_progen != 0, !is.na(second_progen), 
                                                      tstep >= burn_past))
              
              I_coords %>%
                dplyr::select(from = progen_ID, to = ID) %>%
                filter(from != 0) -> edge_df
              gcheck <- graph_from_data_frame(edge_df, directed = TRUE)
              chains <- as.data.frame(list(case_ID = as.numeric(names(components(gcheck)$membership)), 
                                           chain_ID = components(gcheck)$membership))
              I_coords %>%
                left_join(chains, by = c("ID" = "case_ID")) %>%
                mutate(chain_ID = if_else(is.na(chain_ID), 0, chain_ID), 
                       order_id = 1:n(), 
                       step_id = paste0(chain_ID, "_", progen_ID)) %>%
                group_by(chain_ID) %>%
                mutate(step_id = as.numeric(as.factor(step_id))) %>%
                filter(tstep >= burn_past, chain_ID != 0) %>% ## only those past the burn in and vacc periods
                summarize(chain_size = max(step_id)) -> chains
              
              perc_df <- as.data.frame(list(perc = perc_val[i], cov = pup_vacc[j], sim = k,
                                            chain_length = max(chains$chain_size), 
                                            secondaries = secondaries))
            }
})
write.csv(perc_df, "output/perc/perc_vill_flatCov.csv")

### Then just close it out at the end
print("Done:)")
closeCluster(cl) ## for doMPI
mpi.quit()