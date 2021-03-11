# Candidate models (to benchmark) ---------

# sub_cmd:=-t 1 -n 21 -jn fit -wt 1m -sp analysis/scripts/fit.R -md 'gdal' -ar '1-2'

ind <- commandArgs(trailingOnly = TRUE)

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

# set up args
mod_ind <- as.numeric(ind[1])
cand <- fread(fp("analysis/out/fit/candidates.csv"))[mod_ind, ]
nsims <- 1e4

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

# Set up priors
priors <- list(R0 = function(n) exp(rnorm(n, mean = 0, sd = 0.2)),
               iota = function(n) exp(rnorm(n, mean = 0, sd = 1)),
               k = function(n) exp(rnorm(n, mean = 0, sd = 2)))

# Get the startup
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

start_up <- setup_sim(tmax = cand$steps * (cand$nyears),
                      rast = out$rast,
                      death_rate_annual = out$death_rate_annual,
                      birth_rate_annual = out$birth_rate_annual,
                      waning_rate_annual = 1/3,
                      block_fun = block_cells,
                      params = list(start_pop = out$start_pop),
                      step = cand$steps,
                      by_admin = cand$by_admin)

# a little helper function to get names
get_name <- function(i) {
  
  seq <- ifelse(i$sequential, "seq", "kern")
  scale <- ifelse(i$by_admin, "vill", "grid1x1")
  est_incs <- ifelse(i$estincs, "estincs", "fixedincs")
  move <- ifelse(i$weights, "moveprob", "movecont")
  paste(scale, move, est_incs, seq, i$seed, i$partition, sep = "_")
  
}

if(!cand$estincs) {
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
}

if(cand$weights) {
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
disp_fn <- ifelse(cand$sequential, simrabid::steps_weibull, simrabid::dispersal_lognorm)
inc_fn <- ifelse(cand$estincs, simrabid::sim_incursions_pois, simrabid::sim_incursions_hardwired)
move_fn <- ifelse(cand$weights, simrabid::sim_movement_prob, simrabid::sim_movement_continuous)

R0_vals <- priors$R0(nsims)
k_vals <- priors$k(nsims)
if(cand$estincs) iota_vals <- priors$iota(nsims) else iota_vals <- rep(0, nsims)

foreach(j = seq_len(nsims), .combine = rbind, 
        .packages = c("simrabid", "dplyr", 
                      "data.table", "sf", "raster", 
                      "magrittr"), 
        .options.RNG = cand$seed) %dorng% {
     
      simstats <- simrabid(start_up, 
                           start_vacc = cand$start_vacc, 
                           I_seeds = cand$I_seeds, 
                           vacc_dt = vacc_dt,
                           params = c(R0 = R0_vals[j], k = k_vals[j], 
                                      iota = iota_vals[j], 
                                      param_defaults),
                           days_in_step = cand$days_in_step,
                           observe_fun = beta_detect_monthly,
                           serial_fun = serial_lognorm,
                           dispersal_fun = disp_fn, 
                           secondary_fun = nbinom_constrained, 
                           incursion_fun = inc_fn, 
                           movement_fun = move_fn,
                           sequential = cand$sequential, 
                           allow_invalid = cand$allow_invalid,
                           leave_bounds = cand$leave_bounds, 
                           max_tries = 100,
                           summary_fun = test_sim, # use inc_stats after testing!
                           track = cand$track,
                           weights = weights, 
                           row_probs = NULL,
                           coverage = FALSE,
                           break_threshold = cand$break_threshold,
                           by_admin = cand$by_admin)
      
      data.table(simstats, R0 = R0_vals[j], k = k_vals[j], 
                 iota = iota_vals[j])
  } -> out

file_out <- paste0("analysis/out/fit/", get_name(cand), ".csv")

write_create(out,
             fp(file_out),
             fwrite)


# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/fit"

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)
close_cl(cl)

print("Done:)")