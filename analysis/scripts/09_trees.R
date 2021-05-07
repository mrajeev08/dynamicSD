# Get tree stats from posterior sims of best model ----

# sub_cmd:=-t 8 -n 21 -jn tree -wt 1m -md \'gdal\' -ar \'1-2\' -cmd \'1000\'

arg <- commandArgs(trailingOnly = TRUE)

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = TRUE)

# set up args and cluster if applicable ----
nsims <- ifelse(!set_up$slurm, 2, as.numeric(arg[2]))

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
incursions <- fread(fp("analysis/out/incursions.csv"))
load("data/sd_case_data.rda")

# source other scriptss
source("R/sd_data.R")
source("R/get_observed_data.R")
source("R/utils-data.R")
source("R/summ_stats.R")
source("R/run_mods.R")
source("R/score_preds.R")
source("R/conn_metrics.R")

# candidate models to run: the two top choices at the village and grid scale
mod_comp <- readRDS(fp("analysis/out/mod_comp/full_all.rds"))
mod_comp$modlookup %>%
  mutate(modindex = paste0("g", modindex)) %>%
  separate(modname, into = c("scale", "limits", "intros", "move"), sep = "_", 
           remove = FALSE) %>%
  left_join(mod_comp$mod_predicted) %>%
  mutate(N = sum(votes.V1), prob = votes.V1/N) -> lookup

lookup %>%
  group_by(scale) %>%
  arrange(desc(votes.V1)) %>%
  slice(1) -> chosen

int_ind <- ifelse(!set_up$slurm, 2, as.numeric(arg[1]))
if(int_ind == 1) {
  cand_now <- filter(chosen, scale == "grid1x1")$modname
} else {
  cand_now <- filter(chosen, scale == "vill")$modname
}
cands_all <- fread(fp("analysis/out/candidates.csv"))
cands_all$name <- get_name(cands_all, root = TRUE)
cand <- cands_all[name == cand_now][1, ]

# set up best sims ---
# get the startup space
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

post_list <- list.files(fp("analysis/out/par_ests/"))
post_full <- post_list[grep(cand_now, post_list)]
post_full <- post_full[grep("full", post_full)]
post_full <- readRDS(fp(paste0("analysis/out/par_ests/", post_full)))$preds

vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = cand$start_date,
                       date_fun = lubridate::dmy, days_in_step = cand$days_in_step,
                       rollup = 4)

# get observed data ----
sd_case_data %>%
  mutate(cal_month = get_cal_month(dmy(symptoms_started), 
                                   origin_date = cand$start_date)) %$%
  tabulate(cal_month, get_cal_month(cand$apprx_end_date, 
                                    origin_date = cand$start_date)) -> cases_by_month
obs_data <- list(cases_by_month = cases_by_month)

# run simulations ----

set.seed(cand$seed) # same posteriors by candidate
posts <- post_full[, sample(resp, nsims, prob = V1, replace = TRUE), by = "param"]
posts <- split(posts$V1, posts$param)

out_post_sims <- 
  run_simrabid(cand = cand, 
               mod_specs = out,
               param_ests = posts,
               param_defaults = param_defaults,
               nsims = nsims, 
               extra_pars = list(obs_data = obs_data),
               vacc_dt = vacc_dt,
               combine_fun = 'rbind', 
               summary_fun = chain_stats, 
               merge_fun = data.table, 
               secondary_fun = nbinom_constrained,
               weight_covars = list(0), 
               weight_params = list(0), 
               multi = FALSE, 
               convolve_steps = TRUE) 


out_post_sims$modname <- cand_now

# Write out results & close ----
file_pref <- "analysis/out/trees/"

write_create(out_post_sims,
             fp(paste0(file_pref, paste0(cand_now, ".csv"))),
             data.table::fwrite)

# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/trees"

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)
close_cl(cl)

print("Done:)")