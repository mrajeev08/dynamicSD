run_candidate_mod <- function(cand, 
                              priors,
                              nsims, 
                              combine_fun, 
                              name_fun, 
                              seed = 2341) {
  
  
  if(cand$steps == 52) days_in <- 7 else days_in <- 1
  if(sequential) disp <- steps_weibull else disp <- dispersal_lognorm
  
  inc_fun <- sim_incursions_pois 
  
  if(is.na(cand$iota)) {
    param_list <- c(param_list, 
                    list(cell_ids = sample(start_up$cell_ids, tmax, replace = TRUE), 
                         tstep = sample(tmax, tmax, replace = TRUE)))
    
    inc_fun <- sim_incursions_hardwired
  }
  
  if(cand$weights) {
    move <- sim_movement_prob
     } else {
    move <- sim_movement_continuous
    weights <- NULL
  }
  
  foreach(i = seq_len(nsims), .options.RNG = seed) %dorng% {
    
    test <-
      tryCatch(
        expr = {
          simrabid(start_up, 
                   start_vacc = cand$start_vacc, 
                   I_seeds = cand$I_seeds, 
                   vacc_dt = cand$vacc_dt,
                   params = param_list,
                   days_in_step = days_in,
                   observe_fun = beta_detect_monthly,
                   serial_fun = serial_lognorm,
                   dispersal_fun = disp, 
                   secondary_fun = nbinom_constrained, 
                   incursion_fun = inc_fun, 
                   movement_fun = move, 
                   sequential = cand$sequential, 
                   allow_invalid = cand$allow_invalid,
                   leave_bounds = cand$leave_bounds, 
                   max_tries = 100,
                   summary_fun = summary_fun,
                   track = cand$track,
                   weights = weights, 
                   row_probs = NULL,
                   coverage = FALSE,
                   break_threshold = cand$break_threshold,
                   by_admin = admin) 
        },
        error = function(e){
          # write out the error & the parameter log
          err <- cbind(i, e)
          write.csv(err, paste0("logs/err", format(Sys.time(), "%Y%m%d"), ".csv"))
          NULL
        },
        warning = function(w){
          warn <- cbind(i, w)
          write.csv(warn, paste0("logs/warn", format(Sys.time(), "%Y%m%d"), ".csv"))
          NULL
        }
      )
  }
}

