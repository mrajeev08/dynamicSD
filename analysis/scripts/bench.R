# Candidate models (to benchmark) ---------

# sub_cmd:=-t 1 -n 25 -jn times -wt 2m -md 'gdal'

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

# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
load("data/incursions.rda")

# source other scriptss
source("R/sd_data.R")
source("R/utils-data.R")
source("R/summ_stats.R")

# Testing function ---------
vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = "01-Jan-2002",
                       date_fun = lubridate::dmy, units = "weeks", rollup = 4)
# Baseline parameters
pars <- data.frame(track = FALSE,
                   start_vacc = 0.2,
                   break_threshold = 0.85,
                   I_seeds = 0, 
                   death_rate = 0.48, 
                   nyears = 2020 - 2002, 
                   steps = 52, 
                   days_in_step = 7, 
                   leave_bounds = TRUE)

# Set up priors
priors <- list(R0 = function(n) exp(rnorm(n, mean = 0, sd = 0.2)),
               iota = function(n) exp(rnorm(n, mean = 0, sd = 1)),
               k = function(n) exp(rnorm(n, mean = 0, sd = 2)))

# Get the startup
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = 0.48)

start_up_grid <- setup_sim(tmax = pars$steps * (pars$nyears),
                            rast = out$rast,
                            death_rate_annual = out$death_rate_annual,
                            birth_rate_annual = out$birth_rate_annual,
                            waning_rate_annual = 1/3,
                            block_fun = block_cells,
                            params = list(start_pop = out$start_pop),
                            step = pars$steps,
                            by_admin = FALSE)

start_up_vill <- setup_sim(tmax = pars$steps * (pars$nyears),
                           rast = out$rast,
                           death_rate_annual = out$death_rate_annual,
                           birth_rate_annual = out$birth_rate_annual,
                           waning_rate_annual = 1/3,
                           block_fun = block_cells,
                           params = list(start_pop = out$start_pop),
                           step = pars$steps,
                           by_admin = FALSE)

cand <- tidyr::expand_grid(sequential = c(TRUE, FALSE),
                           by_admin = c(TRUE, FALSE),
                           weights = c(TRUE, FALSE), # if prob get weights for each cell!
                           estincs = c(TRUE, FALSE))
cand %<>% 
  mutate(allow_invalid = case_when(weights == TRUE ~ FALSE, 
                                   weights == FALSE ~ TRUE))
cand <- cbind(cand, pars)

# add cell ids & timesteps of incursion to parameter defaults
incursions %<>%
  mutate(cell_id = raster::cellFromXY(out$rast,
                                      cbind(x_coord, y_coord)), 
         tstep = get_timestep(date, 
                              origin_date = "2002-01-01",
                              date_fun = lubridate::ymd, 
                              units = "weeks"))

param_defaults <- c(param_defaults, list(cell_ids = incursions$cell_id, 
                                         tstep = incursions$tstep))

out_times <- foreach(i = iter(cand, by = "row"), .combine = rbind) %do% {
  
  if(i$by_admin) {
    start_up <- start_up_vill
  } else {
    start_up <- start_up_grid
  }
  
  # Weights if using
  if(i$weights) {
    weights <- cell_weights(covars = list(0),
                            params = list(0),
                            start_up$ncell,
                            leave_bounds = TRUE,
                            allow_invalid = FALSE,
                            cells_block = start_up$cells_block,
                            cells_out_bounds = start_up$cells_out_bounds)
  } else {
    weights <- NULL
  }
  
  # params (do this parallel, spit out i + systemtime + wrap in trycatch!
  disp_fn <- ifelse(i$sequential, simrabid::steps_weibull, simrabid::dispersal_lognorm)
  inc_fn <- ifelse(i$estincs, simrabid::sim_incursions_pois, simrabid::sim_incursions_hardwired)
  move_fn <- ifelse(i$weights, simrabid::sim_movement_prob, simrabid::sim_movement_continuous)
  
  nsims <- 1000
  
  R0_vals <- priors$R0(nsims)
  k_vals <- priors$k(nsims)
  iota <- ifelse(i$estincs, priors$iota(nsims), NA)
  
  foreach(j = seq_len(nsims), .combine = 'rbind', 
          .packages = c("simrabid", "dplyr", "data.table", "sf", "raster", 
                        "magrittr")) %dopar% {
                          
    simtm <- system.time(simrabid(start_up, 
                                  start_vacc = i$start_vacc, 
                                  I_seeds = i$I_seeds, 
                                  vacc_dt = vacc_dt,
                                  params = c(R0 = R0_vals[j], k = k_vals[j], 
                                             iota = iota_vals[j], 
                                             param_defaults),
                                  days_in_step = i$days_in_step,
                                  observe_fun = beta_detect_monthly,
                                  serial_fun = serial_lognorm,
                                  dispersal_fun = disp_fn, 
                                  secondary_fun = nbinom_constrained, 
                                  incursion_fun = inc_fn, 
                                  movement_fun = move_fn,
                                  sequential = i$sequential, 
                                  allow_invalid = i$allow_invalid,
                                  leave_bounds = i$leave_bounds, 
                                  max_tries = 100,
                                  summary_fun = test_sim, # use inc_stats after testing!
                                  track = i$track,
                                  weights = weights, 
                                  row_probs = NULL,
                                  coverage = FALSE,
                                  break_threshold = i$break_threshold,
                                  by_admin = i$by_admin))
    c(i, time_est = simtm["elapsed"], R0 = R0_vals[j], k = k_vals[j], 
      iota = iota_vals[j])
  }
  
}

# Output results -----
write_create(out_times,
             fp("analysis/out/test/benchmarks.csv"),
             write.csv, row.names = FALSE)

# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/test"

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)
close_cl(cl)

print("Done:)")