# Comparing Serengeti vax at village vs. district level

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
library(foreach)
library(iterators)
library(doRNG)
library(lubridate)
library(ggplot2)

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
source("R/summ_stats.R")
source("R/run_mods.R")
source("R/conn_metrics.R")

# Testing function ---------
cand <- fread(fp("analysis/out/candidates.csv"))[19, ]
vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = cand$start_date,
                       date_fun = lubridate::dmy, days_in_step = cand$days_in_step,
                       rollup = 4)

# get the startup space
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

start_up <- setup_sim(start_date = cand$start_date, 
                      apprx_end_date =  cand$apprx_end_date, 
                      days_in_step = cand$days_in_step, 
                      rast = out$rast,
                      death_rate_annual = out$death_rate,
                      birth_rate_annual = out$birth_rate,
                      waning_rate_annual = cand$waning,
                      block_fun = block_cells,
                      params = list(start_pop = out$start_pop),
                      by_admin = cand$by_admin)

# removing args if they're included in the priors
pars <- c(R0 = 1.1, k = 1, iota = 1) # just co compare vax!
set.seed(cand$seed)
formals(use_mget)$names <- c("V_mat", "N_mat")

# Randomly allocating
system.time({
  rand <- simrabid(start_up, 
                     start_vacc = cand$start_vacc, 
                     I_seeds = cand$I_seeds, 
                     vacc_dt = vacc_dt,
                     params = c(pars, 
                                param_defaults),
                     days_in_step = cand$days_in_step,
                     observe_fun = beta_detect_monthly,
                     serial_fun = serial_lognorm,
                     dispersal_fun = steps_weibull, 
                     secondary_fun = nbinom_constrained, # function argument 
                     incursion_fun = sim_incursions_pois, 
                     movement_fun = sim_movement_continuous,
                     sequential = TRUE,
                     allow_invalid = TRUE,
                     leave_bounds = TRUE, 
                     max_tries = 100,
                     summary_fun = use_mget, # function argument
                     track = cand$track,
                     weights = NULL, 
                     row_probs = NULL,
                     coverage = FALSE, 
                     break_threshold = cand$break_threshold,
                     by_admin = FALSE)
})
cov_mat_rand <- rand$V_mat/rand$N_mat

# Get district coverage for counterfactual simulation ----
vacc_data <- vacc_dt[, .(vacc_est = sum(vacc_est)), by = "vacc_times"]

vacc_data %>%
  arrange(vacc_times) %>%
  mutate(window = vacc_times - dplyr::lag(vacc_times, 1),
         group = if_else(window <= 8 & !is.na(window), 0, 1),
         group = cumsum(group)) %>%
  group_by(group) %>%
  mutate(vacc_times = min(vacc_times)) %>%
  group_by(vacc_times) %>%
  summarize(vacc_est = sum(vacc_est)) -> vacc_data

pop <- colSums(rand$N_mat)
vacc_data$pop <- pop[vacc_data$vacc_times]
vacc_data$cov <- vacc_data$vacc_est/vacc_data$pop  

vacc_data %>%
  select(vacc_times_dist = vacc_times, cov) %>%
  expand_grid(vacc_locs = 1:75) %>%
  as.data.table(.) -> vacc_data
vacc_data <- vacc_data[vacc_dt, on = "vacc_locs", allow.cartesian = TRUE]
vacc_data %>%
  mutate(vacc_times_diff = abs(vacc_times_dist - vacc_times)) %>%
  group_by(vacc_locs, vacc_times) %>%
  filter(vacc_times_diff == min(vacc_times_diff)) %>%
  select(vacc_locs, vacc_times, vacc_est = cov) %>% 
  mutate(vacc_est = case_when(vacc_times > 200 & vacc_times < 250 ~ vacc_est + 0.2,
                              vacc_times > 400 & vacc_times < 500 ~ vacc_est + 0.17, 
                              vacc_times > 800 & vacc_times < 900 ~ vacc_est + 0.07,
                              TRUE ~ vacc_est)) %>%
  as.data.table(.) -> dist_vacc

fwrite(dist_vacc, "analysis/out/dist_vacc_smoothed.csv")

system.time({
  dist <- simrabid(start_up, 
                   start_vacc = cand$start_vacc, 
                   I_seeds = cand$I_seeds, 
                   vacc_dt = dist_vacc,
                   params = c(pars, 
                              param_defaults),
                   days_in_step = cand$days_in_step,
                   observe_fun = beta_detect_monthly,
                   serial_fun = serial_lognorm,
                   dispersal_fun = steps_weibull, 
                   secondary_fun = nbinom_constrained, # function argument 
                   incursion_fun = sim_incursions_pois, 
                   movement_fun = sim_movement_continuous,
                   sequential = TRUE,
                   allow_invalid = TRUE,
                   leave_bounds = TRUE, 
                   max_tries = 100,
                   summary_fun = use_mget, # function argument
                   track = cand$track,
                   weights = NULL, 
                   row_probs = NULL,
                   coverage = TRUE, 
                   break_threshold = cand$break_threshold,
                   by_admin = FALSE)
})

cov_mat_dist <- dist$V_mat/dist$N_mat

plot(colSums(dist$V_mat)/colSums(dist$N_mat), type = "l", ylim = c(0, 1))
lines(colSums(rand$V_mat)/colSums(rand$N_mat), col = "blue")

hist(cov_mat_dist - cov_mat_rand)
cov_mat_dist <- data.table(cell_id = start_up$cell_ids, cov_mat_dist)
rast_dt <- data.table(as.data.frame(out$rast, xy = TRUE))
rast_dt$cell_id <- 1:nrow(rast_dt)
cov_dt_dist <- rast_dt[cov_mat_dist, on = "cell_id"][, type := "dist"]

cov_mat_rand <- data.table(cell_id = start_up$cell_ids, cov_mat_rand)
cov_dt_rand <- rast_dt[cov_mat_rand, on = "cell_id"][, type := "rand"]

cov_comp <- rbind(cov_dt_dist, cov_dt_rand)

ggplot(cov_comp) + 
  geom_tile(aes(x = x, y = y, fill = V200)) + 
  facet_wrap(~type)

comps_cell <- melt(cov_comp, id.vars = c("type", "x", "y", "layer", "cell_id"))
comps_corr <- dcast(comps_cell, x + y + layer + cell_id + variable ~ type, value.var = "value")

ggplot(comps_corr[!is.na(dist) & dist !=0 & !is.na(rand) & rand != 0]) + 
  geom_hex(aes(dist, rand))
hist(comps_corr$dist - comps_corr$rand)

