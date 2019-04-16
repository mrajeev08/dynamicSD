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
                 vacc_sim = "sim_cell", cov_val = 0.5,
                 return_coords = TRUE)
coords <- check[["I_coords"]]
library(igraph)
coords %>%
  dplyr::select(from = progen_ID, to = ID) %>%
  filter(from != 0) -> edge_df
gcheck <- graph_from_data_frame(edge_df, directed = TRUE)
hist(components(gcheck)$csize)
chains <- as.data.frame(list(case_ID = as.numeric(names(components(gcheck)$membership)), 
                             chain_ID = components(gcheck)$membership))

coords %>%
  left_join(chains, by = c("ID" = "case_ID")) %>%
  mutate(chain_ID = if_else(is.na(chain_ID), 0, chain_ID), 
         order_id = 1:n()) -> check

check %>%
  mutate(observed = rbinom(nrow(.), size = 1, prob = 0.25)) %>%
  filter(observed == 1) %>%
  mutate(progen_obs = ifelse(progen_ID %in% ID, progen_ID, NA)) -> coords_obs50

chain_mat <- matrix(NA, nrow = nrow(coords_obs50), ncol = max(components(gcheck)$membership))
chain_mat[, 1] <- coords_obs50$progen_ID
for(i in 2:ncol(chain_mat)) {
  chain_mat[, i] <- check$progen_ID[match(chain_mat[, i - 1], check$ID)]
}

for(i in 1:nrow(chain_mat)) {
  check$n_unobs[i] <- min(chain_mat[i, ] %in% coords_obs50$ID)
}

n_unobs <- apply(chain_mat, 1, function(x) min(which(x %in% coords_obs50$ID)))
obs_ancestor <- apply(chain_mat, 1, function(x) max(x[x %in% coords_obs50$ID]))

## Next up = algorithm to estimate # of unobserved ancestors and links between cases
## At each step -- assign progen, if prob < x then do two, if prob < x, then do three, 
## if > N progens unobserved call it an incursion?

## Then pipeline to look at across reporting thresholds x vacc thresholds x cut-off thresholds for links & incs

## Need to think about how to define an incursion (closest observed case was an incursion when the ids are -Inf
## or the case itself was an incursion)
