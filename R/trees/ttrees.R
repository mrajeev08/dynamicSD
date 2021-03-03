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
  mutate(observed = rbinom(nrow(.), size = 1, prob = 0.1)) %>%
  filter(observed == 1) %>%
  mutate(progen_obs = ifelse(progen_ID %in% ID, progen_ID, NA)) -> coords_obs50

check %>%
  group_by(chain_ID) %>%
  summarize(chain_size = n()) -> chain_size

chain_mat <- matrix(NA, nrow = nrow(coords_obs50), ncol = max(components(gcheck)$membership))
chain_mat[, 1] <- coords_obs50$progen_ID
for(i in 2:ncol(chain_mat)) {
  chain_mat[, i] <- check$progen_ID[match(chain_mat[, i - 1], check$ID)]
}
coords_obs50$chain_loc <- apply(chain_mat, 1, function(x) length(x[!is.na(x)])) - 1
coords_obs50$n_unobs <- apply(chain_mat, 1, function(x) min(which(x %in% coords_obs50$ID)))
coords_obs50$obs_ancestor <- apply(chain_mat, 1, function(x) max(x[x %in% coords_obs50$ID]))

coords_obs50 %>%
  left_join(chain_size) %>%
  mutate(true_unobs = ifelse(n_unobs == Inf, chain_loc, n_unobs - 1), 
         n_unobs = ifelse(n_unobs == Inf, 0, n_unobs - 1)) -> coords_obs50

nrow(coords_obs50)/(sum(coords_obs50$n_unobs) + nrow(coords_obs50))
nrow(coords_obs50)/(sum(coords_obs50$true_unobs) + nrow(coords_obs50))

library(geosphere)
dist_mat <- as.matrix(dist(cbind(coords_obs50$x_coord, coords_obs50$y_coord)))
t_mat <- outer(coords_obs50$tstep, coords_obs50$tstep, '-')*7

get.SIgamma <- function(nprogens, SIshape = , SIscale = ) {
  
}

get.Kgamma <- function(nprogens, kern_shape = , kern_scale = ) {
  
}

## Options for guess progen function (need a way where one won't be chosen)
## 1. Random uniform guess
## 2. Binomial draw and then random uniform guess (wont always pick one and doesn't make much sense either)
## 3. Multinomial draw based on raw probabilities (will always pick one!)

## 4. Do it as an array, sum across all possible progens?

guess.progen <- function(dist_mat, time_mat, SI_df, kern_df, nmax = 10) {

  for (i in 1:nmax) {
    si_mat <- dgamma(time_mat, SI_df$shape[i], 
                     SI_df$scale[i])
    di_mat <- ifelse(dist_mat <= 100, pgamma(100, kern_df$shape[i], kern_df$scale[i])/100, 
                     dgamma(dist_mat, kern_df$shape[i], kern_df$scale[i]))
    dens_mat <- si_mat*di_mat 
    prob_mat <- dens_mat/rowSums(dens_mat, na.rm = TRUE) 
    progen_mat <- matrix(rbinom(n = nrow(prob_mat)*ncol(prob_mat), size = 1, prob = prob_mat), 
                         nrow = nrow(prob_mat))
    rand <- ravr = runif(1, 0, 1)
    progen_prob <- ifelse(progen_mat == 1, prob_mat, 0)
    progen_prob[is.na(progen_prob)] <- 0
    cumprobs <- cumsum(progen_prob)
    progen_id <- which(cumprobs > ravr, arr.ind = TRUE)[, 1] ## column id of progen 
    logLL <- log(prob_mat[, progen_id])
    found <- rowSums(progen_mat, na.rm = TRUE)
    n_unobs[found != 0, ] <- i
    time_mat[found != 0, ] <- NA ## don't recalculate for ones that were assigned a progenitor
    dist_mat[found != 0, ] <- NA ## don't recalculate for ones that were assigned a progenitor
  }
  
  return(list(n_unobs = n_unobs, progen_id = progen_id, logLL = logLL))
  ## need to create matrices to store logLL at each progen step
}

guess.progen <- function(dist_mat, time_mat, SI_df, kern_df, nmax = 10) {
  n_unobs <- rep(NA, nrow(dist_mat))
  
  for (i in 1:nmax) {
    si_mat <- dgamma(time_mat, SI_df$shape[i], 
                       SI_df$scale[i])
    di_mat <- ifelse(dist_mat <= 100, pgamma(100, kern_df$shape[i], kern_df$scale[i])/100, 
                     dgamma(dist_mat, kern_df$shape[i], kern_df$scale[i]))
    dens_mat <- si_mat*di_mat                
    prob_mat <- dens_mat/rowSums(dens_mat, na.rm = TRUE)
    progen_mat <- matrix(rbinom(n = nrow(prob_mat)*ncol(prob_mat), size = 1, prob = prob_mat), 
                         nrow = nrow(prob_mat))
    rand <- ravr = runif(1, 0, 1)
    progen_prob <- ifelse(progen_mat == 1, prob_mat, 0)
    progen_prob[is.na(progen_prob)] <- 0
    cumprobs <- cumsum(progen_prob)
    progen_id <- which(cumprobs > ravr, arr.ind = TRUE)[, 1] ## column id of progen 
    logLL <- log(prob_mat[, progen_id])
    found <- rowSums(progen_mat, na.rm = TRUE)
    n_unobs[found != 0, ] <- i
    time_mat[found != 0, ] <- NA ## don't recalculate for ones that were assigned a progenitor
    dist_mat[found != 0, ] <- NA ## don't recalculate for ones that were assigned a progenitor
  }
  return(list(n_unobs = n_unobs, progen_id = progen_id, logLL = logLL))
  ## need to create matrices to store logLL at each progen step
}

## Then pipeline to look at across reporting thresholds x vacc thresholds x cut-off thresholds for links & incs

## Need to think about how to define an incursion (closest observed case was an incursion when the ids are -Inf
## or the case itself was an incursion)
