# Simulate from posteriors ---------

# sub_cmd:=-t 2 -n 21 -jn psim -wt 1m -md \'gdal\' -ar \'1-24\' -cmd \'1e3\'

arg <- commandArgs(trailingOnly = TRUE)

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = TRUE)

# set up args and cluster if applicable ----
mod_ind <- ifelse(!set_up$slurm, 1, as.numeric(arg[1]))
nsims <- ifelse(!set_up$slurm, 5, as.numeric(arg[2]))

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
source("R/score_preds.R")

# candidate model to run ----
cands_all <- fread(fp("analysis/out/candidates.csv"))
cands_all$name <- get_name(cands_all, root = TRUE)
cand_now <- unique(cands_all$name)[mod_ind]
# just take the first partition row as these are duplicated
cand <- cands_all[name == cand_now][1, ] 

# get the startup space
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

# loop for independent & joint posteriors ----
post_list <- list.files(fp("analysis/out/par_ests/"))

post_full <- post_list[grep(cand_now, post_list)]
post_full <- post_full[grep("full", post_full)]
post_full <- readRDS(fp(paste0("analysis/out/par_ests/", post_full)))$preds
inds_joint <- post_full[, .(check = all(V1 > 0)), by = "sim_id"][check == TRUE]
post_joint <-  post_full[sim_id %in% inds_joint$sim_id]
post_loop <- list(full = post_full, joint = post_joint)

# loop for vaccination scenarios ----
vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = cand$start_date,
                       date_fun = lubridate::dmy, days_in_step = cand$days_in_step,
                       rollup = 4)
endemic <- vacc_dt[0]
  
sd_vill_pop <- data.table(pop = out$start_pop, 
                          vacc_locs = out$rast[])[, .(pop = sum(pop, na.rm = TRUE)),
                                                  by = "vacc_locs"][!is.na(vacc_locs)]

# loop list
vacc_loop <- list(vill_vacc = vacc_dt, endemic = endemic)

# get observed data ----
sd_case_data %>%
  mutate(cal_month = get_cal_month(dmy(symptoms_started), 
                                   origin_date = cand$start_date)) %$%
  tabulate(cal_month, get_cal_month(cand$apprx_end_date, 
                                    origin_date = cand$start_date)) -> cases_by_month
obs_data <- list(cases_by_month = cases_by_month)

# run simulations ----
out_post_sims <- 
  foreach(i = seq_len(length(vacc_loop)),
          .combine = comb, .multicombine = TRUE) %:%
    foreach(j = seq_len(length(post_loop)), 
            .combine = comb, .multicombine = TRUE) %do% {
      
      # posteriors
      posts <- post_loop[[j]]
      set.seed(cand$seed) # same posteriors by candidate
      posts <- posts[, sample(resp, nsims, prob = V1, replace = TRUE), by = "param"]
      posts <- split(posts$V1, posts$param)
      
      outs <- run_simrabid(cand = cand, 
                           mod_specs = out,
                           param_ests = posts,
                           param_defaults = param_defaults,
                           nsims = nsims, 
                           extra_pars = list(obs_data = obs_data),
                           vacc_dt = vacc_loop[[i]],
                           combine_fun = 'rbind', 
                           summary_fun = ts_stats, 
                           merge_fun = data.table, 
                           secondary_fun = nbinom_constrained,
                           weight_covars = list(0), 
                           weight_params = list(0), 
                           multi = FALSE) 
      
      out_scores <- get_tempstats(outs, obs_data, quants = c(0.5, 0.9), 
                                  nsamp = 100, ncurves = 100, nbest = 5)
      curve_scores <- out_scores$curve_scores
      out_scores <- out_scores$ts_scores
      outs$modname <- out_scores$modname <- curve_scores$modname <- cand_now
      outs$vacc_type <- out_scores$vacc_type <- curve_scores$vacc_type <- names(vacc_loop)[i]
      outs$post_type <- out_scores$post_type <- curve_scores$vacc_type <- names(post_loop)[j]
      
      list(sims = outs, scores = out_scores, curve_scores = curve_scores)
    }
 
# Write out results & close ----

file_pref <- paste0("analysis/out/post_sims/")

write_create(out_post_sims$sims,
             fp(paste0(file_pref, cand_now, "_sims.csv")),
             data.table::fwrite)

write_create(out_post_sims$scores,
             fp(paste0(file_pref, cand_now, "_ts_scores.csv")),
             data.table::fwrite)

write_create(out_post_sims$curve_scores,
             fp(paste0(file_pref, cand_now, "_curve_scores.csv")),
             data.table::fwrite)

# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/post_sims"

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)
close_cl(cl)

print("Done:)")

