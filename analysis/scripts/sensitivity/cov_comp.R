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
load(fp("analysis/out/incursions.csv"))
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

sd_census_data %>%
  mutate(cell_id = cellFromXY(out$rast, xy = cbind(utm_easting, utm_northing)),
         dogs_total = adult_dogs + pups, dogs_vax = adult_dogs_vacc + pups_vacc) %>%
  group_by(cell_id) %>%
  summarize(cov_cell = sum(dogs_vax)/sum(dogs_total)) -> cell_cov
row_probs <- cell_cov$cov_cell[match(start_up$cell_ids, cell_cov$cell_id)]
row_probs[is.na(row_probs) | row_probs == 0] <- 1e-3

system.time({
  prob <- simrabid(start_up, 
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
                     row_probs = row_probs,
                     coverage = FALSE, 
                     break_threshold = cand$break_threshold,
                     by_admin = FALSE)
})

cov_mat_cens <- prob$V_mat/prob$N_mat

hist(cov_mat_cens - cov_mat_rand)
cov_mat_cens <- data.table(cell_id = start_up$cell_ids, cov_mat_cens)
rast_dt <- data.table(as.data.frame(out$rast, xy = TRUE))
rast_dt$cell_id <- 1:nrow(rast_dt)
cov_dt_cens <- rast_dt[cov_mat_cens, on = "cell_id"][, type := "cens"]

cov_mat_rand <- data.table(cell_id = start_up$cell_ids, cov_mat_rand)
cov_dt_rand <- rast_dt[cov_mat_rand, on = "cell_id"][, type := "rand"]

cov_comp <- rbind(cov_dt_cens, cov_dt_rand)

ggplot(cov_comp) + 
  geom_tile(aes(x = x, y = y, fill = V600)) + 
  facet_wrap(~type)

comps_cell <- melt(cov_comp, id.vars = c("type", "x", "y", "layer", "cell_id"))
comps_corr <- dcast(comps_cell, x + y + layer + cell_id + variable ~ type, value.var = "value")

ggplot(comps_corr[!is.na(cens) & cens !=0 & !is.na(rand) & rand != 0]) + 
  geom_hex(aes(cens, rand))
hist(comps_corr$cens - comps_corr$rand)

plot(colSums(rand$V_mat)/colSums(rand$N_mat), type = "l")
lines(colSums(prob$V_mat)/colSums(prob$N_mat), col = "blue")
