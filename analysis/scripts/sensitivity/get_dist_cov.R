# Comparing Serengeti vax by census based probs vs. random allocation

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = TRUE)

cl <- make_cl(set_up$ncores)
register_cl(cl)
print(paste("Cluster size:", cl_size(cl)))

if(!set_up$slurm) fp <- here::here else fp <- cl_safe

# Dependencies
library(raster)
library(data.table)
library(sf)
library(tidyr)
library(dplyr)
library(magrittr)
library(simrabid) # remotes::install_github("mrajeev08/simrabid")

# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
incursions <- fread(fp("analysis/out/incursions.csv"))
load("data/sd_case_data.rda")

# source other scriptss
source("R/sd_data.R")
source("R/get_observed_data.R")
source("R/utils-data.R")

# Testing function ---------
cand <- fread(fp("analysis/out/candidates.csv"))[19, ]
vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = cand$start_date,
                       date_fun = lubridate::dmy, days_in_step = cand$days_in_step,
                       rollup = 4)

vacc_data <- vacc_dt[, .(vacc_est = sum(vacc_est)), by = "vacc_times"]

vacc_data %>%
  arrange(vacc_times) %>%
  mutate(window = vacc_times - dplyr::lag(vacc_times, 1),
         group = if_else(window <= 4 & !is.na(window), 0, 1),
         group = cumsum(group)) %>%
  group_by(group) %>%
  mutate(vacc_times = min(vacc_times)) %>%
  group_by(vacc_times) %>%
  summarize(vacc_est = sum(vacc_est)) -> vacc_data
