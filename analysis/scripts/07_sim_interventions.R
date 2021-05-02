# Simulate from vaccination campaigns -------

# sub_cmd:=-t 6 -n 12 -jn ints -wt 1m -md \'gdal\' -ar \'1-3\' -cmd \'100\' -sn

arg <- commandArgs(trailingOnly = TRUE)

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = FALSE)

# set up args and cluster if applicable ----
vacc_ind <- ifelse(!set_up$slurm, 1, as.numeric(arg[1]))
nsims <- ifelse(!set_up$slurm, 5, as.numeric(arg[2]))
sim_vacc <- "fixed"

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
library(scoringRules)
library(igraph)

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
source("R/score_preds.R")
source("R/conn_metrics.R")

# candidate model to run: top choice from model selection ----
mod_comp <- readRDS(fp("analysis/out/mod_comp/full_grid.rds"))
lookup <- mod_comp$modlookup
mod_comp <- mod_comp$lda_obs$modindex
cand_now <- lookup$modname[match(gsub("g", "", mod_comp), lookup$modindex)][1]

cands_all <- fread(fp("analysis/out/candidates.csv"))
cands_all$name <- get_name(cands_all, root = TRUE)
cand <- cands_all[name == cand_now][1, ]
cand$apprx_end_date <- "2017-12-31" # approximately 15 years 
cand$start_vacc <- 0

# get the startup space
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

# sim from joint posteriors ----
post_list <- list.files(fp("analysis/out/par_ests/"))
post_full <- post_list[grep(cand_now, post_list)]
post_full <- post_full[grep("full", post_full)]
post_full <- readRDS(fp(paste0("analysis/out/par_ests/", post_full)))$preds
inds_joint <- post_full[, .(check = all(V1 > 0)), by = "sim_id"][check == TRUE]
post_joint <-  post_full[sim_id %in% inds_joint$sim_id]

# sim campaigns across range of scenarios ----
# vaccination loop
vacc_scenarios <- expand.grid(vacc_cov = seq(0.2, 1, by = 0.2),
                              vacc_prop = seq(0.2, 1, by = 0.2))
vacc_scenarios <- rbind(data.frame(vacc_cov = 0, vacc_prop = 0), vacc_scenarios)
vacc_scenarios$burn_in <- 5
vacc_scenarios$years <- 10

# 5 years burn in + 5 years vacc + 5 years monitoring
extra_pars <- list(tmin = 10 * 52, rast = out$rast, 
                   adj_dt = get_nb_dt(sd_shapefile))

# setting up other params
int_ind <- ifelse(!set_up$slurm, 1, as.numeric(arg[1]))

# reducing superspreading (as district coverage increased limit transmission)
if(int_ind == 1) {
  nbinom_constrained <- function (n, params = list(R0 = 1.2, 
                                                   k = 1,
                                                   max_secondaries = 100), 
                                  names = c("S", "N")) {
    
    # get S/N 
    list2env(use_mget(names, envir_num = 2), envir = environment())
    
    # prop reduction
    prop_sus <- sum(S)/sum(N)
    if(prop_sus < 0.2) prop_sus <- 0.2
    lims <- params$max_secondaries * prop_sus^2
    secondaries <- rnbinom(n, size = params$k, mu = params$R0)
    secondaries[secondaries > lims] <- round(lims) # out integer 
    return(secondaries)
    
  }
  
  int_name <- "limit_sspreaders" 
} 

# reducing incursions
if(int_ind == 2) {
  sim_incursions_pois <- function(cell_ids, 
                                  params = list(iota = 1),
                                  names = c("S", "N")) {
    
    # get S/N 
    list2env(use_mget(names, envir_num = 2), envir = environment())

    # prop reduction
    prop_sus <- sum(S)/sum(N)
    if(prop_sus < 0.3) prop_sus <- 0.3
    lims <- params$iota * prop_sus^2
    
    n_incs <- rpois(1, lims)
    cell_id <- safe_sample(x = cell_ids, size = n_incs, replace = TRUE)
    
    return(cell_id)
  }
  
  int_name <- "limit_incs" 
  
}

if(int_ind == 3) {
  pup_vacc <- 0.1 # 10 % of pups are vaccinated each week
  int_name <- "limit_temphet" 
} else {
  pup_vacc <- 0
}

out_post_sims <- 
  foreach(i = seq_len(nrow(vacc_scenarios)),
          .combine = rbind) %do% {
            
            # sample from posteriors
            set.seed(cand$seed)
            posts <- post_joint[, sample(resp, nsims, prob = V1, replace = TRUE), by = "param"]
            posts <- split(posts$V1, posts$param)
            
            extra_pars$cov_threshold <- vacc_scenarios[i, "vacc_cov"]
            
            outs <- run_simrabid(cand = cand, 
                                 mod_specs = out,
                                 param_ests = posts,
                                 param_defaults = param_defaults,
                                 nsims = nsims, 
                                 extra_pars = extra_pars,
                                 vacc_dt = vacc_scenarios[i, ],
                                 combine_fun = 'rbind', 
                                 summary_fun = conn_stats, 
                                 merge_fun = data.table, 
                                 secondary_fun = nbinom_constrained,
                                 weight_covars = list(0), 
                                 weight_params = list(0), 
                                 multi = FALSE,
                                 sim_vacc = sim_vacc, 
                                 routine_vacc = pup_vacc) 
            
            outs$cov <- extra_pars$cov_threshold
            outs$prop <- vacc_scenarios[i, "vacc_prop"]
            outs
          }

# Write out results & close ----
file_pref <- paste0("analysis/out/int_sims/", sim_vacc, "_")

write_create(out_post_sims,
             fp(paste0(file_pref, int_name, ".csv")),
             data.table::fwrite)

# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/int_sims"

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)
close_cl(cl)

print("Done:)")