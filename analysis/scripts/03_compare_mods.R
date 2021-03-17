# Candidate models (to benchmark) ---------

# sub_cmd:=-t 1 -n 21 -jn fit -wt 1m -md 'gdal' -ar '1-4' -cmd '100'

arg <- commandArgs(trailingOnly = TRUE)

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

# set up args
mod_ind <- as.numeric(arg[1])
cand <- fread(fp("analysis/out/fit/candidates.csv"))[mod_ind, ]
nsims <- as.numeric(arg[2])

# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
load("data/incursions.rda")
load("data/sd_case_data.rda")

# source other scriptss
source("R/sd_data.R")
source("R/get_observed_data.R")
source("R/utils-data.R")
source("R/summ_stats.R")
source("R/run_mods.R")

# Testing function ---------
vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = cand$start_date,
                       date_fun = lubridate::dmy, days_in_step = cand$days_in_step,
                       rollup = 4)

# get the startup space
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

# Set up priors
priors <- list(R0 = function(n) exp(rnorm(n, mean = 0.1, sd = 0.2)), # centered around 1.1
               iota = function(n) exp(rnorm(n, mean = 0.5, sd = 0.5)), # centered around 1.5
               k = function(n) exp(rnorm(n, mean = 0.5, sd = 0.5))) # uniform

# get observed data
obs_data <- get_observed_data(sd_case_data, 
                              cand = cand, 
                              out = out)$obs_data

out_sims <- run_simrabid(cand = cand, 
                         mod_specs = out,
                         param_ests = priors,
                         param_defaults = param_defaults,
                         nsims = nsims, 
                         extra_pars = list(obs_data = obs_data),
                         vacc_dt = vacc_dt,
                         combine_fun = 'rbind', 
                         summary_fun = inc_stats, 
                         secondary_fun = nbinom_constrained,
                         weight_covars = list(0), 
                         weight_params = list(0))  

file_out <- paste0("analysis/out/fit/", get_name(cand), ".csv")


write_create(out_sims,
             fp(file_out),
             data.table::fwrite)


# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/fit"

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)
close_cl(cl)

print("Done:)")