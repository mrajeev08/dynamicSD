# Candidates --------

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = TRUE)

# Dependencies
library(data.table)
library(tidyr)
library(dplyr)
library(magrittr)

# List of candidate models to output ---------
# Baseline parameters
pars <- data.frame(track = FALSE,
                   start_vacc = 0.2,
                   break_threshold = 0.85,
                   I_seeds = 0, 
                   death_rate = 0.465, 
                   days_in_step = 7, 
                   waning = 1/3,
                   start_date = "2002-01-01",
                   apprx_end_date = "2020-12-31")

# candidate df to write out
cand <- tidyr::expand_grid(sequential = c(TRUE, FALSE),
                           by_admin = c(TRUE, FALSE),
                           weights = c(TRUE, FALSE), # if prob get weights for each cell!
                           estincs = c(TRUE, FALSE), 
                           partition = c(1, 2), 
                           leave_bounds = c(TRUE, FALSE)) # 32 rows
cand %<>% 
  mutate(allow_invalid = case_when(weights == TRUE ~ FALSE, 
                                   weights == FALSE ~ TRUE), 
         exclude = !weights & !leave_bounds) %>%
  filter(!exclude)

cand <- cbind(cand, pars)

set.seed(124)
cand$seed <- sample(1e4, size = nrow(cand), replace = FALSE)

# Output results -----
write_create(cand,
             here::here("analysis/out/candidates.csv"),
             fwrite, row.names = FALSE)

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)

