### Functions for running most pared down IBM possible
## Malavika Rajeev
## December 2018

# Setting up to run in parallel
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
res(r) <- 1000 # 1 x 1km 
SD_raster <- rasterize(SD_shape, r)

## Data by cellID
grid_data <- get.demdata(shapefile = SD_shape, res_meters = 1000, pop = pop_data, census = census_data)
values(SD_raster) <- grid_data$cell_id ## cell IDs as values

# values(SD_raster) <- grid_data$start_pop ## cell IDs as values
# plot(SD_raster, bty = "n")
# myCol <- colorRampPalette(c("pink", "red", "darkred"))
# 
# # Add axis labels
# plot(SD_raster, 
#      col = myCol(6),
#      axes = FALSE, box = FALSE)

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

## function for combining foreach outputs
multicomb <- function(x, ...) {
  mapply(cbind, x, ..., SIMPLIFY = FALSE)
}

##' SD scenario ------------------------------------------------------------------------------------------------
##' scaling factors for incursions
scale_iota <- read_csv("data/d_relative_incursion_rates.csv")
scale_iota <- scale_iota$relative.number
iota = 1; 
scale_iota = c(rep(scale_iota, each = 4), 
               rep(scale_iota[length(scale_iota)], 4));
nsim <- 10
dir_name <- paste0("output/simsSD")

if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

SD_sims <- foreach(i = 1:nsim, .combine = multicomb, .multicombine = TRUE,
                   .packages = c("data.table", "sp", "raster")) %dopar% {
  check <- sim.IBM(R_0 = 1.1, k = 0.6, inc_rate = 1*scale_iota, p_revacc = 1, start_vacc = 0.2, 
                   I_seed = 1, p_vacc = grid_data$cov, return_coords = FALSE)
  N <- check[["N"]]
  S <- check[["S"]]
  I_all <- check[["I_all"]]
  I_ts <- colSums(matrix(colSums(I_all), nrow = 4)) ## monthly ts
  S_ts <- colMeans(matrix(colSums(S), nrow = 4)) 
  N_ts <- colMeans(matrix(colSums(N), nrow = 4)) 
  list(I_ts = I_ts, S_ts = S_ts, N_ts = N_ts)
}

write.csv(SD_sims[["I_ts"]], paste0(dir_name, "/Imat_SDvacc.csv"))
write.csv(SD_sims[["S_ts"]], paste0(dir_name, "/Smat_SDvacc.csv"))
write.csv(SD_sims[["N_ts"]], paste0(dir_name, "/Nmat_SDvacc.csv"))

##' Endemic ------------------------------------------------------------------------------------------------
vacc_mat <- matrix(0, nrow = nrow(grid_data), ncol = tmax)
SD_sims <- foreach(i = 1:nsim, .combine = multicomb, .multicombine = TRUE,
                   .packages = c("data.table", "sp", "raster")
                   ) %dopar% {
  check <- sim.IBM(R_0 = 1.1, k = 0.5, inc_rate = 1*scale_iota, p_revacc = 1, start_vacc = 0, 
                   p_vacc = grid_data$cov, return_coords = FALSE)
  N <- check[["N"]]
  S <- check[["S"]]
  I_all <- check[["I_all"]]
  I_ts <- colSums(matrix(colSums(I_all), nrow = 4)) ## monthly ts
  S_ts <- colMeans(matrix(colSums(S), nrow = 4)) 
  N_ts <- colMeans(matrix(colSums(N), nrow = 4)) 
  list(I_ts = I_ts, S_ts = S_ts, N_ts = N_ts)
}

write.csv(SD_sims[["I_ts"]], paste0(dir_name, "/Imat_SDendemic.csv"))
write.csv(SD_sims[["S_ts"]], paste0(dir_name, "/Smat_SDendemic.csv"))
write.csv(SD_sims[["N_ts"]], paste0(dir_name, "/Nmat_SDendemic.csv"))


##' Elimination ----------------------------------------------------------------------------------------
##' ## Need for time step functions within data cleaning
burn_in <- 5
sim_yrs <- 15
tmax <- (burn_in + sim_yrs) * 52
burn_past <- tmax - burn_in*52

SD_sims <- foreach(i = 1:nsim, .combine = multicomb, .multicombine = TRUE,
                   .packages = c("igraph", "data.table", "dplyr", "sp", "raster")) %dopar% {
  vacc_mat <- sim.campaigns(data = grid_data, vills = unique(grid_data$villcode),
                            sim_years = sim_yrs, burn_in_years = burn_in, 
                            vill_weeks = sample(1:52, 75, replace = TRUE))
  
  check <- sim.IBM(R_0 = 1.1, k = 0.5, inc_rate = 1, p_revacc = 1, I_seed = 1,
                   p_vacc = grid_data$cov, start_vacc = 0, perc_val = 1,
                   vacc_sim = "sim_vill", cov_val = 0.7,
                   return_coords = FALSE)
  N <- check[["N"]] 
  S <- check[["S"]]
  I_all <- check[["I_all"]]
  I_ts <- colSums(matrix(colSums(I_all), nrow = 4)) ## monthly ts
  S_ts <- colMeans(matrix(colSums(S), nrow = 4)) 
  N_ts <- colMeans(matrix(colSums(N), nrow = 4)) 
  list(I_ts = I_ts, S_ts = S_ts, N_ts = N_ts)
}

write.csv(SD_sims[["I_ts"]], paste0(dir_name, "/Imat_SDelim.csv"))
write.csv(SD_sims[["S_ts"]], paste0(dir_name, "/Smat_SDelim.csv"))
write.csv(SD_sims[["N_ts"]], paste0(dir_name, "/Nmat_SDelim.csv"))

### Then just close it out at the end
print("Done:)")
closeCluster(cl) ## for doMPI
mpi.quit()

##' Looking at cov distribution
##' ----------------------------------------------------------------------------------------
##' ## Need for time step functions within data cleaning
burn_in <- 0
sim_yrs <- 15
tmax <- (burn_in + sim_yrs) * 52
burn_past <- tmax - burn_in*52
vacc_mat <- sim.campaigns(data = grid_data, vills = unique(grid_data$villcode),
                                               sim_years = sim_yrs, burn_in_years = burn_in, 
                                               vill_weeks = sample(1:52, 75, replace = TRUE))
                     
check <- sim.IBM(R_0 = 1.1, k = 0.5, inc_rate = 1, p_revacc = 1, I_seed = 1,
                  p_vacc = grid_data$cov, start_vacc = 0, perc_val = 1,
                  vacc_sim = "sim_vill", cov_val = 0.7,
                  return_coords = FALSE)
N <- check[["N"]] 
S <- check[["S"]]
cov_mat <- 1 - S/N
sims_SD <- data.table::melt(as.data.table(cov_mat), variable.name = "loc", value.name = "cov")
ggplot(data = sims_SD, aes(x = cov)) +
  geom_histogram(binwidth = 0.1, color = "white", fill = "grey") +
  xlab("Coverage") +
  ylab("Frequency")
